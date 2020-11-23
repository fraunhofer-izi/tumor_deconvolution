#! /usr/bin/env Rscript
suppressPackageStartupMessages(library(ArrayExpress))
suppressPackageStartupMessages(library(SummarizedExperiment))

outPath = "../data/rnaseq-R"
outPath = normalizePath(outPath)

assay(experimentSummary$rnaseq)

sets = queryAE(keywords = "endothelia", species = "homo+sapiens")
AEset = ArrayExpress("E-MEXP-1416")
