#! /usr/bin/env Rscript
suppressPackageStartupMessages({
    library(utils)
    library(foreach)
    library(doMC)
})
registerDoMC(32)

with_array = TRUE
rescale_arrays = FALSE
rescale_array_to_n_reads = 5e7
important_genes_table = "tabula-SupplementaryTables.csv"
importance = 2 # weight factor for important genes
gene_remove = c(coarse = .8, fine = .8)
workPath = "../data"
workPath = normalizePath(workPath)
inRDS = file.path(workPath, "rnaseq-R-HUGO-GSE")
inRDSarray = file.path(workPath, "rnaseq-R-HUGO-GSE-arrays")
outName = paste0("tissue_tsv_HUGO_imp", ifelse(with_array, "_array", ""))
outTSV = file.path(workPath, outName)
HUGOfile = file.path(workPath, "HUGO_mpa.RDS")
fl_name = file.path(inRDS, "info.RDS")
array_fl_name = file.path(inRDSarray, "info.RDS")
dir.create(outTSV, showWarnings = FALSE)

# read data
HUGO_map = readRDS(HUGOfile)
imp_genes = read.csv(important_genes_table, comment.char="#")

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
features = unique(HUGO_map$HUGO)

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
stat_cols = function(info_df, set_name, cell_type) {
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

comb = function(...) mapply(cbind, ..., SIMPLIFY = FALSE)

# getting gene statistics
gstats = foreach (grain = names(info)) %do% {
    message(paste("Unifying", grain, "..."))
    cell_types = setdiff(names(info[[grain]]),
                         c("NA", NA, "NA.", "too.unspecific"))
    stats = foreach(cell_type = cell_types, .combine = cbind) %:%
        when(cell_type != "others") %dopar% {
        dat_list = info[[grain]][[cell_type]]
        ind = sapply(dat_list, is.null)
        cstats = foreach(set_name = names(dat_list[!ind])) %dopar%
            stat_cols(dat_list[[set_name]], set_name, cell_type)
        res = do.call(comb, cstats)
        good_feature = (res[["hits"]] - res[["min_missing"]]) >= 1
        good_feature = good_feature & (res[["max_ambiguity"]] == 1)
        NAs = good_feature * res[["NAs"]]
        total_good_feature = good_feature * res[["total"]]
        total_samps = sum(res[["total"]][1, ])
        quality = apply(total_good_feature - NAs, 1, sum) / total_samps
        qdf = data.frame(quality)
        colnames(qdf) = cell_type
        return(qdf)
    }
    min_quality = apply(stats, 1, min)
    iind = names(min_quality) %in% important_genes
    min_quality[iind] = min_quality[iind] * importance
    threshold = quantile(min_quality, gene_remove[grain])
    ind = min_quality >= threshold
    message(paste("Keeping", sum(ind), "features..."))
    message(paste("...and", sum(iind&ind), "of the", length(important_genes),
                  "mapped genes from 'The Cancer Immunome Atlas'."))
    sfeatures = names(min_quality)[ind]
    foreach(cell_type = cell_types) %dopar% {
        dat_list = info[[grain]][[cell_type]]
        ind = sapply(dat_list, is.null)
        if (all(ind)) {
            expr = data.frame(row.names = sfeatures)
            total = 0
        } else {
            expr = foreach(set_name = names(dat_list[!ind]),
                           .combine = cbind) %dopar%
                filtered_expr(set_name, grain, cell_type, sfeatures)
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
