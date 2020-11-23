#! /usr/bin/env Rscript
library(utils)
library(foreach)
library(doMC)
library(ggplot2)
library(bigmemory)

source("R/utilities.R")
source("R/grains.R")
registerDoMC(32)

# Produces feature statistics about extracted data from the GEO studies.

gene_share = c(coarse = .8, fine = .2)

CRC_files = paste0("/<censored_path>/dominik.otto/",
                   "gdc-data/colorectal_cancer/")
BRCA_files = paste0("/<censored_path>/dominik.otto/",
                    "gdc-data/breast-cancer/")


inPath = "../data/rnaseq-R"
outPath = "../data/"
inPath = normalizePath(inPath)
outPath = normalizePath(outPath)
inFile = file.path(outPath, "extraction_stats.RDS")
outRDS = file.path(outPath, "rnaseq-R")
outTSV = file.path(outPath, "tissue_tsv")

shared = big.matrix(nrow = 1, ncol = 1, type = 'integer')
task = quote({
    expr = NULL
    if (dir.exists(fl)) {
        htseq_file = dir(fl, full.names = TRUE)
        htseq_file = grep(".htseq.counts.*gz$", htseq_file, value = TRUE)
        if (length(htseq_file) != 1) stop(paste("Ambigious:", fl))
        zfile = gzfile(htseq_file, "rt")
        expr = read.table(zfile, header = FALSE, row.names = 1)
        close(zfile)
        colnames(expr) = basename(fl)
    }
    shared[1] <- shared[1] + 1
    setTxtProgressBar(pb, shared[1])
    return(expr)
})

message('Collecting CRC experessions ...')
CRC_cases = dir(CRC_files, full.names = TRUE)
pb = txtProgressBar(min = 1, max = length(CRC_cases), style = 3)
shared[1] = 0
CRC_features = foreach(fl = CRC_cases,
                       .combine = merge_all) %dopar% eval(task)
close(pb)
saveRDS(CRC_features, file.path(outRDS, "all_CRC.RDS"))
no_version = gsub("\\.[^.]+$", "", rownames(CRC_features))
rownames(CRC_features) = no_version
meta = grepl("^_", no_version)
CRC_features = subset(CRC_features, !meta)
write.table(CRC_features, file.path(outTSV, "CRC.tsv"),
            quote = FALSE, fileEncoding = "utf-8", sep = "\t")

message('Collecting BRCA experessions ...')
BRCA_cases = dir(BRCA_files, full.names = TRUE)
pb = txtProgressBar(min = 1, max = length(BRCA_cases), style = 3)
shared = big.matrix(nrow = 1, ncol = 1, type = 'integer')
shared[1] = 0
BRCA_features = foreach(fl = BRCA_cases,
                        .combine = merge_all) %dopar% eval(task)
close(pb)
saveRDS(BRCA_features, file.path(outRDS, "all_BRCA.RDS"))
no_version = gsub("\\.[^.]+$", "", rownames(BRCA_features))
rownames(BRCA_features) = no_version
meta = grepl("^_", no_version)
BRCA_features = subset(BRCA_features, !meta)
write.table(BRCA_features, file.path(outTSV, "BRCA.tsv"),
            quote = FALSE, fileEncoding = "utf-8", sep = "\t")


message('Collecting cell type experessions ...')
selected_features = rownames(CRC_features)
chall_tab = readRDS(inFile)
chall_tab = subset(chall_tab, ENS_genes & raw_counts > 0 & !strange)

for (grain in c('coarse', 'fine')) {
    message(paste("Doing", grain, "grain ..."))
    grainf = paste0(grain, '.cell.type')
    cell.types = setdiff(unique(grains[, grain]),
                         "too unspecific")
    cell.types = make.names(cell.types)

    write_tables = function(selected_features, rm.na = FALSE) {
        missings = foreach(ct = cell.types) %dopar% {
            ct_tab = chall_tab[chall_tab[, grainf] == ct, ]
            expr = data.frame(row.names = selected_features)
            if (nrow(ct_tab) != 0) {
                GSEs = unique(ct_tab$gse)
                expr = foreach(GSE = GSEs, .combine = cbind) %do% {
                    selection = subset(ct_tab, gse == GSE)
                    fl_name = paste0(GSE, "_raw-counts.RDS")
                    expr_path = file.path(inPath, fl_name)
                    expr = tryCatch(readRDS(expr_path),
                        error = function(e) return(NULL))
                    if (is.null(expr)) return(NULL)
                    no_subsamp = gsub("\\.[^.]+$", "", colnames(expr))
                    ind = no_subsamp %in% selection$id
                    no_version = gsub("\\.[^.]+$", "", rownames(expr))
                    map = match(selected_features, no_version)
                    expr = expr[map, ind]
                    rownames(expr) = selected_features
                    return(expr)
                }
                is_integer = is.int(expr, na.rm = TRUE)
                totals = apply(expr, 2, sum, na.rm = TRUE)
                positive = apply(expr>0, 2, all, na.rm = TRUE)
                not_normalized = (positive & totals > 15e5) | is_integer
                expr = expr[, not_normalized]
                # Since names are made unique per GSE this removes
                # only samples that occure in multiple GSEs:
                expr = expr[, !duplicated(colnames(expr))]
            }
            fl_name = file.path(outRDS, paste0(grain, "-",
                                               make.names(ct), ".RDS"))
            if(isFALSE(rm.na)) saveRDS(expr, fl_name)
            expr = expr[, !duplicated(colnames(expr))]
            if (isTRUE(rm.na)) {
                na_sample = apply(is.na(expr), 2, any)
                expr = expr[, !na_sample]
                missings = table(na_sample)
            } else {
                missings = apply(is.na(expr), 1, sum)
                missings["total"] = ncol(expr)
            }
            fl_name = file.path(outTSV, paste0(grain, "-",
                                               make.names(ct), ".tsv"))
            write.table(expr, fl_name, quote = FALSE,
                        fileEncoding = "utf-8", sep = "\t")
            message(paste("Saved", ct, "with dimensions",
                          paste(dim(expr), collapse="-")))
            return(missings)
        }
        names(missings) = cell.types
        return(missings)
    }

    missings = write_tables(selected_features)
    names(missings) = NULL
    mv = do.call(c, missings)
    mstat = aggregate(mv, list(gene = names(mv)), FUN = sum)
    bound = quantile(mstat$x, gene_share[grain])
    sf = mstat$gene[mstat$x <= bound]

    # save reduced sets
    message("Removing all NAs ...")
    missing_types = write_tables(sf, rm.na = TRUE)
    crc_file = paste0(grain, "-CRC.tsv")
    write.table(CRC_features[sf, ], file.path(outTSV, crc_file),
                quote = FALSE, fileEncoding = "utf-8", sep = "\t")
    brca_file = paste0(grain, "-brca.tsv")
    write.table(BRCA_features[sf, ], file.path(outTSV, brca_file),
                quote = FALSE, fileEncoding = "utf-8", sep = "\t")

}
