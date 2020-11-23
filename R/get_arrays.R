#! /usr/bin/env Rscript
message("Loading packages ...")
suppressPackageStartupMessages({
    library(GEOquery)
    library(GEOmetadb)
    library(AnnotationDbi)
    library(R.utils)
    library(data.table)
    library(bigmemory)
    library(foreach)
    library(doMC)
})
registerDoMC(32)
source("R/grains.R")

grain_levels = c("coarse", "fine")
outPath = "../data/"
outRDS = file.path(outPath, "rnaseq-R-HUGO-GSE-arrays")
dir.create(outRDS, showWarnings = FALSE)
HUGOfile = file.path(outPath, "HUGO_mpa.RDS")
info_fl_name = file.path(outRDS, "info.RDS")
chall_tab = read.csv2("geo-expression-array-immune-cells.csv")
HUGO_map = readRDS(HUGOfile)
gene_cols = c("Gene Symbol", "gene_assignment", "GENE_SYMBOL",
              "GENE_ID", "UCSC_RefGene_Name", "ID")

for (grain in grain_levels) {
    map = match(chall_tab$cell.type, grains$original)
    grainf = paste0(grain, '.cell.type')
    chall_tab[, grainf] = grains[map, grain]
}

gses = unique(chall_tab$gse)
GEOmetadbFileName = "GEOmetadb.sqlite"
geoSQLfile = file.path(outPath, GEOmetadbFileName)
if (!file.exists(geoSQLfile)) {
    message("Getting GEOmetadb SQL file ....")
    geoSQLfile = getSQLiteFile(destdir = outPath)
}
sqlcon = dbConnect(SQLite(), geoSQLfile)

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

message("Downloading and mapping micro arrays ...")
shared = big.matrix(nrow = 1, ncol = 1, type = 'integer')
shared[1] = 0
pb = txtProgressBar(min = 1, max = length(gses), style = 3)

info = foreach(gse = gses, .combine = info_combine) %dopar% {
    capture.output({capture.output({
        dat <- getGEO(gse, GSEMatrix = TRUE)
    }, file = "/dev/null", type = "message")
    }, file = "/dev/null", type = "output")
    info = foreach(p = seq_along(dat), .combine = info_combine) %do% {
        d = dat[[p]]
        expr = exprs(d)
        sub_tab = subset(chall_tab, id %in% colnames(expr))
        if (nrow(sub_tab) == 0) return(NULL)
        fd = fData(d)
        if (is.null(ncol(fd)) || ncol(fd) == 0) return(NULL)
        probes = as.character(fd[, 1])
        gpl = annotation(d)
        query = paste0("select bioc_package from gpl where gpl='",
                       gpl, "'")
        pck = NA
        try({pck <- dbGetQuery(sqlcon, query)$bioc_package})
        if (is.na(pck) || pck == "NA" || length(probes) == 0)
            return(NULL)
        if (!grepl(".db$", pck)) pck = paste0(pck, ".db")
        suppressPackageStartupMessages({
            if (!require(pck, character.only = TRUE, quietly = TRUE)) {
                capture.output({capture.output({
                    renv::install(pck)
                }, file = "/dev/null", type = "message")
                }, file = "/dev/null", type = "output")
            }
        })
        anno = eval(parse(text = paste0(pck, "::", pck)))
        if (!"SYMBOL" %in% columns(anno)) return(NULL)
        if (!any(probes %in% keys(anno))) return(NULL)
        suppressMessages(genes <- mapIds(anno, k = probes,
                                         column = c("SYMBOL"),
                                         keytype = "PROBEID"))
        map = match(toupper(genes), HUGO_map$ref, nomatch = 1)
        mapped_features = HUGO_map[map, ]
        rownames(mapped_features) = rownames(expr)

        info = foreach(grain = grain_levels) %do% {
            grainf = paste0(grain, '.cell.type')
            cell.types = make.names(unique(sub_tab[, grainf]))
            info = foreach(cell.type = cell.types) %do% {
                gsms = sub_tab$id[make.names(sub_tab[, grainf]) == cell.type]
                samps = as.character(gsms)
                gset = expr[, samps, drop = FALSE]
                dt_expr = data.table(gset)[, lapply(.SD, as.numeric)]
                na_dt = dt_expr[, lapply(.SD, is.na)]
                na_dt$HUGO = mapped_features$HUGO
                na_dt = na_dt[, lapply(.SD, sum), by="HUGO"]
                dt_expr = cbind(dt_expr, mapped_features)
                dt_expr = dt_expr[, c(list(hits = .N,
                                           max_ambiguity = max(matches)),
                                   lapply(.SD, mean, na.rm = TRUE)),
                                   by = "HUGO", .SDcols = samps]
                info_df = data.frame(
                    dt_expr[, c("HUGO", "hits", "max_ambiguity")])
                fexpr = data.frame(
                    dt_expr[, ..samps], row.names = dt_expr$HUGO)
                info_df$min_missing = apply(na_dt[, ..samps], 1, min)
                info_df$NAs = apply(is.na(fexpr), 1, sum)
                info_df$total = ncol(fexpr)

                gse_p = paste0(gse, "-p", p)
                fl_name = paste0(outRDS, "/", gse_p, "-", grain,
                                 "-", cell.type, ".RDS")
                saveRDS(fexpr, fl_name)
                info = list()
                info[[gse_p]] = info_df
                return(info)
            }
            names(info) = cell.types
            return(info)
        }
        names(info) = grain_levels
        return(info)
    }
    shared[1] = shared[1] + 1
    setTxtProgressBar(pb, shared[1])
    return(info)
}
close(pb)
dbDisconnect(sqlcon)
message("Saving mapping info ...")
saveRDS(info, info_fl_name)

