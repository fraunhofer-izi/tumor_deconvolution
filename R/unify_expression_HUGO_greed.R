#! /usr/bin/env Rscript
#SBATCH --job-name="optimize feature selection"
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=10-00:00:00
#SBATCH --output=logs/feature_selection.%j.slurmlog

suppressPackageStartupMessages({
    library(utils)
    library(foreach)
    library(doMC)
    library(bigmemory)
})
options(bigmemory.typecast.warning = FALSE)
registerDoMC(16)

with_array = TRUE
rescale_arrays = FALSE
rescale_array_to_n_reads = 5e7
quality_threshold = 5
important_genes_table = "tabula-SupplementaryTables.csv"
start_features_file = "start_features.txt"
importance = 2 # weight factor for important genes
gene_remove = c(coarse = .8, fine = .8)
workPath = "../data"
workPath = normalizePath(workPath)
inRDS = file.path(workPath, "rnaseq-R-HUGO-GSE")
inRDSarray = file.path(workPath, "rnaseq-R-HUGO-GSE-arrays")
outName = paste0("tissue_tsv_greed", ifelse(with_array, "_array", ""))
outTSV = file.path(workPath, outName)
HUGOfile = file.path(workPath, "HUGO_mpa.RDS")
fl_name = file.path(inRDS, "info.RDS")
array_fl_name = file.path(inRDSarray, "info.RDS")
annoFile = file.path(workPath, "sample_data.tsv")
dir.create(outTSV, showWarnings = FALSE)

# read data
HUGO_map = readRDS(HUGOfile)
imp_genes = read.csv(important_genes_table, comment.char="#")
anno = read.delim(annoFile)
start_features = readLines(start_features_file)

message("Loading and combining mapping info ...")
info = readRDS(fl_name)
info_combine = function(A, B) {
    if (!is.list(A) || !is.list(B)) return(c(A, B))
    names(A) = make.names(names(A))
    names(B) = make.names(names(B))
    info = list()
    for (k in union(names(A), names(B))) {
        if (!k %in% names(A)) {
            info[[k]] = B[[k]]
        } else if (!k %in% names(B)) {
            info[[k]] = A[[k]]
        } else {
            info[[k]] = info_combine(A[[k]], B[[k]])
        }
    }
    return(info)
}
info = info_combine(info, readRDS(array_fl_name))
features = setdiff(as.character(unique(HUGO_map$HUGO)), "none")

# getting important genes of The Cancer Immunome Atlas (https://tcia.at/)
imp_genes$Metagene = sub("ClQA", "C1QA", imp_genes$Metagene) # fix typo
imp_genes$Metagene = sub("ClQB", "C1QB", imp_genes$Metagene) # fix typo
map = match(toupper(imp_genes$Metagene), toupper(HUGO_map$ref))
nmap = is.na(map)
if (any(nmap)) {
    genes = paste0(imp_genes[nmap, "Metagene"],
                   " (", imp_genes[nmap, "Cell.type"], ")",
                   collapse="; ")
    warning(paste("The following", sum(nmap),
                 "genes of 'The Cancer Immunome Atlas' could not be mapped to",
                 "HUGO gene names:", genes))
}
important_genes = as.character(HUGO_map$HUGO[na.omit(map)])

# defining helper functions
stat_cols = function(info_df, set_name, cell_type, features = features) {
    total = info_df$total[1]
    map = match(features, info_df$HUGO)
    rdf = data.frame(info_df[map, ], row.names = features)
    ind = is.na(map)
    if (any(ind)) {
        rdf[ind, ] = data.frame(HUGO = features[ind], hits = 0,
                                min_missing = total,
                                max_ambiguity = 0,
                                NAs = total, total = total)
    }
    extract = function(vec) {
        mdf = data.frame(as.vector(vec), row.names = features)
        colnames(mdf) = paste0(set_name, "_", cell_type)
        return(mdf)
    }
    result = lapply(rdf, extract)
    result[["HUGO"]] = NULL
    return(result)
}

normalize_array = function(x) {
    x = x - min(x, na.rm = TRUE)
    if (isTRUE(rescale_arrays)) {
        total = sum(x, na.rm = TRUE)
        x = x * (rescale_array_to_n_reads / total)
    }
    return(x)
}

filtered_expr = function(set_name, grain, cell_type, features) {
    set_name = gsub("\\.p", "-p", set_name)
    is_array = grepl("^[^-]*-p[0-9]*$", set_name)
    inDir = ifelse(is_array, inRDSarray, inRDS)
    fl_name = file.path(inDir, paste0(set_name, "-",
                                       grain, "-",
                                       make.names(cell_type), ".RDS"))
    if (!file.exists(fl_name)) return(NULL)
    expr = readRDS(fl_name)
    expr = expr[features, , drop = FALSE]
    if (isTRUE(is_array)) {
        expr = apply(expr, 2, normalize_array)
    }
    rownames(expr) = features
    return(expr)
}

filtered_samples_expr = function(set_name, grain, cell_type, features) {
    expr = filtered_expr(set_name, grain, cell_type, features)
    local_gse = gsub("\\..*", "", set_name)
    local_anno = subset(anno, gse == local_gse)
    map = match(gsub("\\..*", "", colnames(expr)), local_anno$id)
    good_sample = local_anno[map, "aggregated_quality"] > quality_threshold
    good_sample[is.na(good_sample)] = TRUE
    return(expr[, good_sample, drop = FALSE])
}

score_total = function(elements, features) {
    required_dat = elements[features, , drop = FALSE]
    bad_samples = apply(required_dat == 0, 2, any)
    return(sum(required_dat[, !bad_samples]))
}

get_sample_wise = function(elements, features) {
    required_dat = elements[features, , drop = FALSE]
    bad_samples = apply(required_dat == 0, 2, any)
    return(apply(required_dat[, !bad_samples, drop = FALSE], 2, sum))
}

score_add_total = function(base, new) {
    ind = new > 0
    return(sum(base[ind]) + sum(new[ind]))
}

logsumexp = function (x) {
    y = max(x)
    return(y + log(sum(exp(x - y))))
}

# getting gene statistics
weight = rep(1, length.out = length(features))
weight[features %in% important_genes] = 2
shared = big.matrix(nrow = 1, ncol = 1, type = "integer")
gstats = foreach (grain = names(info)) %do% {
    message(paste("Unifying", grain, "..."))
    cell_types = setdiff(names(info[[grain]]),
                         c("NA", NA, "NA.", "too.unspecific", "others"))
    avail = foreach(cell_type = cell_types) %dopar% {
        dat_list = info[[grain]][[cell_type]]
        ind = sapply(dat_list, is.null)
        if (all(ind)) {
            message(paste("No data for", cell_type))
            return(NULL)
        }
        result = foreach(set_name = names(dat_list[!ind]),
                       .combine = cbind, .inorder = FALSE) %do% {
            expr = filtered_samples_expr(set_name, grain, cell_type, features)
            if (nrow(expr) == 0) {
                return(NULL)
            }
            stats = stat_cols(dat_list[[set_name]], set_name, cell_type, features)
            good_feature = (stats[["hits"]] - stats[["min_missing"]]) >= 1
            good_feature = good_feature & (stats[["max_ambiguity"]] == 1)
            return((!is.na(expr) + 0) * (weight + drop(good_feature)))
        }
        return(result)
    }
    names(avail) = cell_types
    selection = features %in% start_features
    names(selection) = features
    message("Obtimizing feature selection ...")
    avail_shared = list()
    for (cell_type in cell_types) {
        avail_shared[[cell_type]] = as.big.matrix(avail[[cell_type]], type = "integer")
    }
    current_score = NA
    iteration = 1
    report = data.frame()
    while(TRUE) {
        shared[1] = 0
        pb = txtProgressBar(min = 1, max = length(features), style = 3)
        blocks = lapply(avail, get_sample_wise, selection)
        scores = foreach(f = as.character(features), .combine = rbind,
                         .inorder = TRUE) %dopar% {
            local_selection = selection
            if (local_selection[f]) {
                local_selection[f] = FALSE
                type_scores = foreach(av = avail_shared, .combine = c) %do% {
                    return(score_total(av, local_selection))
                }
            } else {
                type_scores = foreach(ct = names(blocks), .combine = c) %do% {
                    sample_names = names(blocks[[ct]])
                    return(score_add_total(blocks[[ct]],
                                           avail_shared[[ct]][f, sample_names]))
                }
            }
            report = data.frame(t(type_scores))
            colnames(report) = cell_types
            report$total = -logsumexp(-type_scores/10)
            report$feature = f
            shared[1] = shared[1] + 1
            setTxtProgressBar(pb, shared[1])
            return(report)
        }
        close(pb)
        next_selection = selection
        next_selection[next_feature] = !next_selection[next_feature]
        nfeatures = sum(next_selection)
        feature_scores = scores$total
        next_feature = which.max(feature_scores)
        next_score = feature_scores[next_feature]
        scores$nfeatures = nfeatures
        report = rbind(report, scores[next_feature, ])
        message(paste(date(), "iteration:", iteration, "- genes:", nfeatures,
                      " - score:", next_score))
        if (!is.na(current_score) && next_score < current_score) {
            message("Converged.")
            break
        }
        selection = next_selection
        current_score = next_score
    }
    iind = names(selection) %in% important_genes
    message(paste("Keeping", sum(selection), "features..."))
    message(paste("...and", sum(iind&selection), "of the", length(important_genes),
                  "mapped genes from 'The Cancer Immunome Atlas'."))
    saveRDS(report, file.path(outTSV,  paste0(grain, "-report.RDS")))
    sfeatures = names(selection[selection])
    foreach(cell_type = cell_types) %dopar% {
        dat_list = info[[grain]][[cell_type]]
        ind = sapply(dat_list, is.null)
        if (all(ind)) {
            expr = data.frame(row.names = sfeatures)
            total = 0
        } else {
            expr = foreach(set_name = names(dat_list[!ind]),
                           .combine = cbind) %dopar%
                filtered_samples_expr(set_name, grain, cell_type, sfeatures)
            total = ncol(expr)
            good_samps = !apply(is.na(expr), 2, any)
            good_samps = good_samps & !duplicated(colnames(expr))
            expr = expr[sfeatures, good_samps, drop = FALSE]
        }
        message(paste("Keeping", ncol(expr), "of", total, cell_type))
        fl_name = file.path(outTSV, paste0(grain, "-",
                                           make.names(cell_type), ".tsv"))
        write.table(expr, fl_name, quote = FALSE,
                    fileEncoding = "utf-8", sep = "\t")
    }
}
