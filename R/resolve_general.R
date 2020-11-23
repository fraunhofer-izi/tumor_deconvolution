#! /usr/bin/env Rscript

outPath = "../data/rnaseq-R"
outPath = normalizePath(outPath)
dir.create(outPath, showWarnings = FALSE)

source("R/resolvers.R")

arg = commandArgs(trailingOnly=TRUE)
inPath = arg[1]
ncores = 1
if (length(arg)>1) ncores = as.integer(arg[2])
GSE_number = basename(inPath)
message(paste("Doing", GSE_number))
if (ncores > 1) registerDoMC(ncores)

resolvers = c("htseq" = htseq_extract,
              "kallisto" = kallisto_extract,
              "featureCounts" = featureCounts_extract,
              "samplemap" = samplemap_extract,
              "mmseq" = mmseq_extract,
              "excel" = excel_extract,
              "ENSGlist" = ENSGlist_extract)

resolved = FALSE
if (GSE_number %in% exclude) {
    expr = NULL
    resolved = TRUE
} else if (GSE_number %in% names(custom_resolvers)) {
    expr = custom_resolvers[[GSE_number]]()
    resolved = TRUE
} else {
    files = dir(inPath, full.names = TRUE)
    if (length(files) == 0) {
        expr = NULL
        resolved = TRUE
    } else {
        formats = get_format(files)
        for (r in names(resolvers)) {
            if (r %in% formats) {
                ind = formats == r
                expr = resolvers[[r]](files[ind])
                resolved = TRUE
                break
            }
        }
    }
}

if (!resolved) {
    if (any(formats == "unrecognized")){
        stop(paste("No resolver for", GSE_number))
    } else {
        message(paste("No extractable format found for", GSE_number))
        expr = NULL
    }
}

if (!is.null(expr)) {
    outFile = paste0(GSE_number, "_raw-counts")
    outFile = file.path(outPath, outFile)
    saveRDS(expr, file = paste0(outFile, ".RDS"))
    write.table(expr, file = paste0(outFile, ".tsv"),
                quote = FALSE, sep = '\t')
    message(paste("Resolved", GSE_number))
}
