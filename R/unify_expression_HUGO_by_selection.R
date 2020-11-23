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

if (grepl("/R$", getwd())) setwd("..")

with_array = TRUE
rescale_arrays = FALSE
rescale_array_to_n_reads = 5e7
quality_threshold = 5
important_genes_table = "tabula-SupplementaryTables.csv"
feature_file = c(coarse = "../data/feature_wsy_selection_coarse.txt",
                 fine = "../data/feature_wsy_selection_fine.txt")
workPath = "../data"
workPath = normalizePath(workPath)
inRDS = file.path(workPath, "rnaseq-R-HUGO-GSE")
inRDSarray = file.path(workPath, "rnaseq-R-HUGO-GSE-arrays")
outName = paste0("tissue_tsv_wsa", ifelse(with_array, "_array", ""))
outTSV = file.path(workPath, outName)
HUGOfile = file.path(workPath, "HUGO_mpa.RDS")
array_fl_name = file.path(inRDSarray, "info.RDS")
rnaseq_fl_name = file.path(inRDS, "info.RDS")
annoFile = file.path(workPath, "sample_data.tsv")
dir.create(outTSV, showWarnings = FALSE)

# read data
HUGO_map = readRDS(HUGOfile)
anno = read.delim(annoFile)
imp_genes = read.csv(important_genes_table, comment.char="#")

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

message("Loading mapping info ...")
info_combine = function(A, B) {
    if (!is.list(A) || !is.list(B)) return(c(A, B))
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
info = info_combine(readRDS(array_fl_name), readRDS(rnaseq_fl_name))

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

gstats = foreach (grain = names(info)) %do% {
    message(paste("Unifying", grain, "..."))
    selection = readLines(feature_file[grain])
    nimp = sum(important_genes %in% selection)
    message(paste("Keeping", length(selection), "features..."))
    message(paste("...and", nimp, "of the", length(important_genes),
                  "mapped genes from 'The Cancer Immunome Atlas'."))
    cell_types = setdiff(names(info[[grain]]),
                         c("NA", NA, "NA.", "too.unspecific", "others"))
    dummy = foreach(cell_type = cell_types) %do% {
        dat_list = info[[grain]][[cell_type]]
        expr = foreach(set_name = names(dat_list),
                       .combine = cbind) %dopar%
            filtered_samples_expr(set_name, grain, cell_type, selection)
        total = ncol(expr)
        good_samps = !apply(is.na(expr), 2, any)
        expr = expr[selection, good_samps, drop = FALSE]
        message(paste("Keeping", ncol(expr), "of", total, cell_type))
        fl_name = file.path(outTSV, paste0(grain, "-",
                                           make.names(cell_type), ".tsv"))
        write.table(expr, fl_name, quote = FALSE,
                    fileEncoding = "utf-8", sep = "\t")
    }
}
