#! /usr/bin/env Rscript
library(utils)
library(foreach)
library(doMC)
library(data.table)
registerDoMC(8)
setDTthreads(1)

source("R/resolvers.R")
source("R/grains.R")

grain_levels = c("coarse", "fine")
inPath = "../data/rnaseq-R"
outPath = "../data/"
inPath = normalizePath(inPath)
outPath = normalizePath(outPath)
#inFile = file.path(outPath, "extraction_stats.RDS")
inFile = file.path(outPath, "sample_data.tsv")
outRDS = file.path(outPath, "rnaseq-R-HUGO-GSE")
HUGOfile = file.path(outPath, "HUGO_mpa.RDS")
temp_wd = file.path(outPath, "rslurm_temp")
info_fl_name = file.path(outRDS, "info.RDS")
dir.create(temp_wd, showWarnings = FALSE)
dir.create(outRDS, showWarnings = FALSE)

added_gse = as.character(read.csv2("geo-additions.csv")$gse)
HUGO_map = readRDS(HUGOfile)
#chall_tab = readRDS(inFile)
chall_tab = read.delim(inFile)

make_RDS = function(set_name, cell_type, ct_tab) {
    fl_name = paste0(set_name, "_raw-counts.RDS")
    expr_path = file.path(inPath, fl_name)
    expr = tryCatch(suppressWarnings(readRDS(expr_path)),
        error = function(e) return(NULL))
    if (is.null(expr)) return(NULL)
    if (set_name %in% names(filter_additionals)) {
        expr = filter_additionals[[set_name]](expr, cell_type)
    } else if (set_name %in% ct_tab$gse) {
        not_unique = gsub("\\..*$", "", colnames(expr))
        ind = not_unique %in% ct_tab$id
        if (!any(ind)) return(NULL)
        expr = expr[, ind, drop = FALSE]
    }
    if (is.null(expr) || ncol(expr) == 0 || nrow(expr) == 0) return(NULL)
    positive = apply(expr >= 0, 2, all, na.rm = TRUE)
    if (!any(positive)) return(NULL)
    expr = expr[, positive, drop = FALSE]
    samps = colnames(expr)
    features = rownames(expr)
    if (any(grepl("^ENS", head(features)))) {
        features = gsub("\\..*$", "", features)
    }
    map = match(toupper(features), HUGO_map$ref, nomatch = 1)
    mapped_features = HUGO_map[map, ]
    rownames(mapped_features) = rownames(expr)
    dt_expr = data.table(expr)[, lapply(.SD, as.numeric)]
    na_dt = dt_expr[, lapply(.SD, is.na)]
    na_dt$HUGO = mapped_features$HUGO
    na_dt = na_dt[, lapply(.SD, sum), by="HUGO"]
    dt_expr = cbind(dt_expr, mapped_features)
    dt_expr = dt_expr[, c(list(hits = .N, max_ambiguity = max(matches)),
                       lapply(.SD, sum, na.rm = TRUE)),
                       by = "HUGO", .SDcols = samps]
    info = data.frame(dt_expr[, c("HUGO", "hits", "max_ambiguity")])
    expr = data.frame(dt_expr[, ..samps], row.names = dt_expr$HUGO)
    info$min_missing = apply(na_dt[, ..samps], 1, min)
    info$NAs = apply(is.na(expr), 1, sum)
    info$total = ncol(expr)
    fl_name = file.path(outRDS, paste0(make.names(set_name), "-",
                                       grain, "-",
                                       make.names(cell_type), ".RDS"))
    saveRDS(expr, fl_name)
    return(info)
}

info = foreach (grain = grain_levels) %do% {
    message(paste("HUGO mapping", grain))
    grainf = paste0(grain, '.cell.type')
    cell.types = setdiff(unique(grains[, grain]),
                         "too unspecific")
    cell.types = make.names(cell.types)
    info = foreach(cell_type = cell.types) %dopar% {
        ct_tab = chall_tab[chall_tab[, grainf] == cell_type, ]
        sets = c(as.character(unique(ct_tab$gse)), added_gse)
        info = foreach(set = sets) %dopar% make_RDS(set, cell_type, ct_tab)
        names(info) = sets
        return(info)
    }
    names(info) = cell.types
    return(info)
}
names(info) = grain_levels
saveRDS(info, info_fl_name)

#source("R/unify_expression_HUGO.R")
