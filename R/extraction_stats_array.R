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
inFile = "geo-expression-array-immune-cells.csv"
inPath = "../data/rnaseq-R-HUGO-GSE-arrays"
outPath = "../data/"
inPath = normalizePath(inPath)
outPath = normalizePath(outPath)
outFile = file.path(outPath, "extraction_stats_array.RDS")

chall_table = read.csv2(inFile)
GSEs = sort(unique(chall_table$gse))
files = dir(inPath)

pb = txtProgressBar(min = 1, max = length(GSEs), style = 3)
shared = big.matrix(nrow = 1, ncol = 1, type = 'integer')
shared[1] = 0

dat = foreach(GSE = GSEs, .combine = rbind) %dopar% {
    res = subset(chall_table, gse == GSE)
    res$resolved = any(grepl(GSE, files))
    shared[1] = shared[1] + 1
    setTxtProgressBar(pb, shared[1])
    return(res)
}
setTxtProgressBar(pb, length(GSEs))
close(pb)

occurences = table(dat$id)
dat$occurences = occurences[dat$id]

nodup = function(index) {
    index[index] = !duplicated(dat$id[index])
    return(index)
}
dat$status = "duplicates"
ind = nodup(dat$resolved)
dat$status[ind] = "resolved"
uind = nodup(!ind)
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

saveRDS(dat, outFile)
write.table(dat, gsub(".RDS$", ".tsv", outFile), quote = FALSE,
            fileEncoding = "utf-8", sep="\t")

pl = ggplot(dat, aes(x = factor(1), fill = state)) +
    theme_void() + theme(axis.text.x=element_blank()) +
    geom_bar(position="fill") + coord_polar(theta = "y") +
    ggtitle("GEO Expression Array Samples")
    #geom_text(aes(label = ..count..), stat = "count",
    #          position = position_fill(vjust = 0.5))
ggsave(file.path(outPath, "extraction_status_array.png"),
       plot = pl, device = "png", width = 10, height = 7)

raw_dat = subset(dat, status == "resolved")
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
ggsave(file.path(outPath, "original_cell_types_array.png"),
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
ggsave(file.path(outPath, "fine_cell_types_array.png"),
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
ggsave(file.path(outPath, "coarse_cell_types_array.png"),
       plot = ccpl, device = "png", width = 9, height = 7)
