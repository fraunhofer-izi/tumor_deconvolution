#! /usr/bin/env Rscript
library(utils)
library(foreach)
library(ggplot2)
library(doMC)
library(bigmemory)

source("R/utilities.R")
source("R/grains.R")

registerDoMC(4)

# Produces statistics about extracted data from the GEO studies.

verbose = FALSE
min_test_samp = 2
min_test_share = .01
inFile = "geo-rnaseq-immune-cells.csv"
inPath = "../data/rnaseq-R"
outPath = "../data/"
inPath = normalizePath(inPath)
outPath = normalizePath(outPath)
outFile = file.path(outPath, "extraction_stats.RDS")

chall_table = read.csv2(inFile)
GSEs = sort(unique(chall_table$gse))

pb = txtProgressBar(min = 1, max = length(GSEs), style = 3)
shared = big.matrix(nrow = 1, ncol = 1, type = 'integer')
shared[1] = 0

dat = foreach(GSE = GSEs, .combine = rbind) %dopar% {
    res = subset(chall_table, gse == GSE)
    res$resolved = FALSE
    res$repeats = 0
    res$integer_counts = 0
    res$raw_counts = 0
    res$ENS_genes = FALSE
    res$ENS_trans = FALSE
    res$strange = FALSE
    fl_name = paste0(GSE, "_raw-counts.RDS")
    expr_path = file.path(inPath, fl_name)
    expr = NULL
    num.expr = NULL
    if (file.exists(expr_path)) {
        tryCatch({
            expr = readRDS(expr_path)
            dim_string = paste(dim(expr), collapse="-")
            mess = paste(GSE, "has dimensions", dim_string)
            if (isTRUE(verbose)) message(mess)
            res$resolved = TRUE
        }, error = function(e) {})
    }
    if (!is.null(expr)) {
        samps = gsub("\\..*", "", colnames(expr))
        avails = res$id[res$available]
        res$repeats = sapply(res$id, n_hits, b = samps)
        if (any(!sapply(expr, is.numeric))) {
            num.expr = foreach(co = expr, .combine = cbind) %do% {
                as.numeric(as.character(co))
            }
        } else {
            num.expr = expr
        }
        is_integer = is.int(num.expr, na.rm = TRUE)
        int_samps = samps[is_integer]
        res$integer_counts = sapply(res$id, n_hits, b = int_samps)
        totals = apply(num.expr, 2, sum, na.rm = TRUE)
        positive = apply(num.expr>0, 2, all, na.rm = TRUE)
        not_normalized = (positive & totals > 15e5) | is_integer
        raw_samps = samps[not_normalized]
        res$raw_counts = sapply(res$id, n_hits, b = raw_samps)
        ens = grepl("^ENSG[0-9]*", rownames(expr))
        res$ENS_genes = sum(ens, na.rm = TRUE) / length(ens) > .2
        ens = grepl("^ENST[0-9]*", rownames(expr))
        res$ENS_trans = sum(ens, na.rm = TRUE) / length(ens) > .2
        res$strange = any(grepl(";[0-9.]*;", rownames(expr)))
    }
    rm(expr, num.expr)
    shared[1] = shared[1] + 1
    setTxtProgressBar(pb, shared[1])
    return(res)
}

occurences = table(dat$id)
dat$occurences = occurences[dat$id]

nodup = function(index) {
    index[index] = !duplicated(dat$id[index])
    return(index)
}
dat$status = "duplicates"
ind = nodup(dat$raw_counts > 0)
dat$status[ind] = "only gene symbols for raw counts"
eind = ind & (dat$ENS_genes | dat$ENS_trans)
dat$status[eind] = "raw counts w/ ensemble IDs"
nind = nodup(dat$repeats > 0 & !ind)
dat$status[nind] = "only normalized counts"
rind = nodup(dat$resolved & (!ind) & (!nind))
dat$status[rind] = "no data available"
uind = nodup((!ind) & (!nind) & (!rind))
dat$status[uind] = "unresolved"
res_stats = table(dat$status)
n = res_stats[dat$status]
total = length(dat$status)
dat$state = paste0(dat$status, " (n=", n, ", ",
                   signif(100*n/total, 2), "%)")

map = match(dat$cell.type, grains$original)
dat$coarse.cell.type = make.names(grains$coarse[map])
dat$fine.cell.type = make.names(grains$fine[map])
dat$combined.cell.type = make.names(grains$combined[map])

raw_ind = dat$status == "raw counts w/ ensemble IDs" & dat$fine.cell.type != "unresolved"
dat$test_data = FALSE
set.seed(1)
for (cell_type in unique(dat$coarse.cell.type)) {
    ct_ind = dat$coarse.cell.type == cell_type
    avail_ind = raw_ind & ct_ind
    min_share = ceiling(min_test_share * sum(avail_ind))
    nsamp = max(min_test_samp, min_share)
    test_set = sample(which(avail_ind), nsamp)
    dat$test_data = dat$test_data | (1:nrow(dat) %in% test_set)
}

saveRDS(dat, outFile)
write.table(dat, gsub(".RDS$", ".tsv", outFile), quote = FALSE,
            fileEncoding = "utf-8", sep="\t")

pl = ggplot(dat, aes(x = factor(1), fill = state)) +
    theme_void() + theme(axis.text.x=element_blank()) +
    geom_bar(position="fill") + coord_polar(theta = "y") +
    ggtitle("GEO RNA-Seq Samples")
    #geom_text(aes(label = ..count..), stat = "count",
    #          position = position_fill(vjust = 0.5))
ggsave(file.path(outPath, "extraction_status.png"),
       plot = pl, device = "png", width = 10, height = 7)

raw_dat = subset(dat, status == "raw counts w/ ensemble IDs")
cell_stats = table(raw_dat$cell.type)
missings = setdiff(dat$cell.type, raw_dat$cell.type)
missings = paste(missings, "(n=0, 0%)")
n = cell_stats[raw_dat$cell.type]
total = length(raw_dat$cell.type)
raw_dat$ct = paste0(raw_dat$cell.type, " (n=", n,
                           ", ", signif(100*n/total, 2), "%)")
limits = sort(union(raw_dat$ct, missings))
cpl = ggplot(raw_dat, aes(x = factor(1), fill = ct)) +
    theme_void() + theme(axis.text.x=element_blank()) +
    geom_bar(position = "fill") + coord_polar(theta = "y") +
    ggtitle("Fine Grain Cell Types for Raw Counts with Ensemble IDs") +
    scale_fill_discrete(name = "origina cell types",limits = limits)
ggsave(file.path(outPath, "original_cell_types.png"),
       plot = cpl, device = "png", width = 14, height = 7)

cell_stats = table(raw_dat$fine.cell.type)
missings = setdiff(make.names(grains$fine), raw_dat$fine.cell.type)
if (length(missings)>0) missings = paste(missings, "(n=0, 0%)")
n = cell_stats[raw_dat$fine.cell.type]
total = length(raw_dat$fine.cell.type)
raw_dat$cct = paste0(raw_dat$fine.cell.type, " (n=", n,
                           ", ", signif(100*n/total, 2), "%)")
limits = sort(union(raw_dat$cct, missings))
ccpl = ggplot(raw_dat, aes(x = factor(1), fill = cct)) +
    theme_void() + theme(axis.text.x=element_blank()) +
    geom_bar(position = "fill") + coord_polar(theta = "y") +
    ggtitle("Fine Grain Cell Types for Raw Counts with Ensemble IDs") +
    scale_fill_discrete(name = "fine cell types", limits = limits)
ggsave(file.path(outPath, "fine_cell_types.png"),
       plot = ccpl, device = "png", width = 9, height = 7)

cell_stats = table(raw_dat$coarse.cell.type)
missings = setdiff(make.names(grains$coarse), raw_dat$coarse.cell.type)
if (length(missings)>0) missings = paste(missings, "n=0")
n = cell_stats[raw_dat$coarse.cell.type]
total = length(raw_dat$coarse.cell.type)
raw_dat$cct = paste0(raw_dat$coarse.cell.type, " (n=", n,
                           ", ", signif(100*n/total, 2), "%)")
limits = sort(union(raw_dat$cct, missings))
ccpl = ggplot(raw_dat, aes(x = factor(1), fill = cct)) +
    theme_void() + theme(axis.text.x=element_blank()) +
    geom_bar(position = "fill") + coord_polar(theta = "y") +
    ggtitle("Coarse Grain Cell Types for Raw Counts with Ensemble IDs") +
    scale_fill_discrete(name = "coarse cell types", limits = limits)
ggsave(file.path(outPath, "coarse_cell_types.png"),
       plot = ccpl, device = "png", width = 9, height = 7)
