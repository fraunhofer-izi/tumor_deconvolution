#! /usr/bin/env Rscript
library(foreach)
library(doMC)
library(ggplot2)

source("R/utilities.R")
registerDoMC(32)

# Produces feature statistics about extracted data from the GEO studies.

inPath = "../data/rnaseq-R"
outPath = "../data/"
inPath = normalizePath(inPath)
outPath = normalizePath(outPath)
inFile = file.path(outPath, "extraction_stats.RDS")
outFile = file.path(outPath, "feature_stats.RDS")

dat = readRDS(inFile)
ENSind = dat$ENS_genes & !dat$strange
GSEs = unique(dat$gse[ENSind])

agg_features = function(a, b) {
    if (is.null(a) && is.null(b)) {
        return(NULL)
    } else if (is.null(a)) {
        return(b)
    } else if (is.null(b)) {
        return(a)
    }
    dat = merge_all(a, b)
    count = apply(dat, 1, sum, na.rm = TRUE)
    res = data.frame(count = count, row.names = rownames(dat))
    return(res)
}

feature_counts = foreach(GSE = GSEs, .combine = agg_features) %dopar% {
    fl_name = paste0(GSE, "_raw-counts.RDS")
    expr_path = file.path(inPath, fl_name)
    if (file.exists(expr_path)) {
        expr = readRDS(expr_path)
    } else {
        return(NULL)
    }
    no_versions = gsub("\\.[^\\.]+$", "", rownames(expr))
    count = apply(!is.na(expr), 1, sum)
    count = aggregate(count, list(feature = no_versions), max)
    res = data.frame(count = count$x, row.names = count$feature)
    return(res)
}

pl = ggplot(feature_counts, aes(x = count)) + geom_histogram()
print(pl)
