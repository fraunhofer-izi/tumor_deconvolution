#! /usr/bin/env Rscript
suppressPackageStartupMessages(library(utils))
suppressPackageStartupMessages(library(GEOquery))
suppressPackageStartupMessages(library(doMC))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(openxlsx))

studies = "../data/rnaseq-raw-studies"
softs = "../data/rnaseq-raw-studies-soft"
studies = normalizePath(studies)
softs = normalizePath(softs)

anno_table = c(hg19 = paste0("/<censored_path>/dominik.otto/",
                             "annotation-tables/results/",
                             "all_gene_annotations.RData"),
               hg38 = paste0("/<censored_path>/dominik.otto/",
                             "annotation-tables/results/",
                             "all_gene_annotations_V10.RData"))

not_a_sample = c("Length", "length")

n_hits = function(a, b) sum(a == b, na.rm = TRUE)

attempt_load = function(file_name) {
    state = "error"
    data_env = new.env()
    tryCatch({
        load(file_name, envir = data_env)
        state = "loaded"
    })
    if (state == "error") {
        tryCatch({
            data_env$dat = readRDS(file_name)
            state = "loaded"
        })
    }
    if (state == "error") {
        stop(paste("Error loading", file_name))
    }
    return(data_env)
}

map_hg19symbols = function(syms){
    load(anno_table["hg19"]) # loads annoTab
    symbol_map = match(syms, annoTab$geneName)
    gencode_cols = grep("gencode", colnames(annoTab),
                        value = TRUE, ignore.case = TRUE)
    latest_gencode = sort(gencode_cols, decreasing = TRUE)[1]
    new_names = annoTab[symbol_map, latest_gencode]
    new_names = remove_ensg_version(new_names)
    return(new_names)
}

remove_ensg_version = function(ensg) {
    return(gsub("\\.[0-9]*$", "", ensg))
}

is.binary <- function(filepath, max = 1000){
    fl = file(filepath, "rb", raw= TRUE)
    bytes = readBin(fl, "int", max, size = 1, signed = FALSE)
    close(fl)
    return(max(bytes) > 128)
}

get_format = function(files) {
    resolve_with = foreach(fl = files, .combine = c) %do% {
        if (grepl("\\.bed$", fl)) return("bed")
        if (grepl("\\.xlsx?$", fl)) return("excel")
        if (grepl("peaks\\.txt$", fl)) return("peaks")
        if (grepl("\\.CGmap$", fl)) return("CGmap")
        if (grepl("\\.gtf$", fl)) return("gtf")
        if (grepl("/filelist.txt$", fl)) return("filelist")
        bn = basename(fl)
        if (grepl("GPL10558", bn)) return("microarray")
        if (is.binary(fl)) {
            if (grepl("\\.bw", fl)) {
                return("bigWig")
            } else {
                return("binary")
            }
        }
        txt = readLines(fl, 10)
        if (grepl("^track", txt[1])) return("track")
        if (txt[1] == "feature\tcount") return("htseq")
        if (txt[1] == "ID\tRawCounts") return("htseq")
        if (grepl("^[^#\t]+\t+[0-9.]+\t*$", txt[1])) return("htseq")
        if (grepl("featureCounts", txt[1])) return("samplemap")
        pat = "^target_id\tlength\teff_length\test_counts\ttpm$"
        if (grepl(pat, txt[1])) return("kallisto")
        pat = paste0("^gene_id\ttranscript_id(s)\tlengthlength\t",
                     "expected_count\tTPM\tFPKM$")
        if (grepl(pat, txt[1])) return("kallisto")
        if (grepl("^# Mapped fragments: ", txt[1])) return("mmseq")
        pat = "^ENSG[0-9]*(\\.[0-9]*)?[ \t][0-9.]*$"
        if (grepl(pat, txt[2]) && grepl("^GSM[0-9]*_", bn)) return("ENSGlist")
        pat = "^[^\\s#]+[\t,;][0-9.]+[\t,;]"
        ind = grepl(pat, paste0(txt, ";"))
        for (tx in txt[ind]) {
            nf = length(strsplit(tx, "\t")[[1]])
            nf = max(nf, length(strsplit(tx, ",")[[1]]))
            nf = max(nf, length(strsplit(tx, ";")[[1]]))
            if (nf > 2) return("samplemap")
        }
        return("unrecognized")
    }
    return(resolve_with)
}

get_soft = function(gse, verbose = FALSE) {
    soft_path = file.path(softs, gse,
                          paste0(gse, "_family.soft"))
    if (isTRUE(verbose)) {
        soft = getGEO(filename = soft_path)
    } else {
        soft = suppressMessages(getGEO(filename = soft_path))
    }
    return(soft)
}

ENSGlist_extract = function(files) {
    if (length(files) == 0) return(NULL)
    pat = "^ENSG[0-9]*(\\.[0-9]*)?[ \t][0-9.]*$"
    expr = foreach(fl = files, .combine = merge_all) %dopar% {
        first_line = readLines(fl, 1)
        skip = ifelse(grepl(pat, first_line), 0, 1)
        read.table(fl, row.names = 1, skip = skip)
    }
    gsms = gsub("_.*", "", basename(files))
    colnames(expr) = gsms
    return(expr)
}

htseq_extract = function(files) {
    if (length(files) == 0) return(NULL)
    headers = c("feature\tcount", "ID\tRawCounts")
    expr = foreach(fl = files, .combine = merge_all) %dopar% {
        first_line = readLines(fl, 1)
        skip = ifelse(first_line %in% headers, 1, 0)
        read.table(fl, row.names = 1, skip = skip)
    }
    gsms = gsub("_.*", "", basename(files))
    colnames(expr) = gsms
    return(expr)
}

featureCounts_extract = function(files) {
    if (length(files) == 0) return(NULL)
    expr = foreach(fl = files, .combine = cbind) %dopar% {
        read.table(fl, row.names = 1, skip = 1, header = T)[6]
    }
    gsms = gsub("_.*", "", files)
    colnames(expr) = gsms
    return(expr)
}

merge_all = function(A, B) {
    samps = c(colnames(A), colnames(B))
    A[, "rn"] = rownames(A)
    B[, "rn"] = rownames(B)
    res = join(A, B, by = "rn", type = "full", match = "first")
    rownames(res) = res$rn
    return(res[, samps, drop = FALSE])
}

samplemap_extract = function(files) {
    if (length(files) == 0) return(NULL)
    GSE = basename(dirname(files[1]))
    soft = get_soft(GSE)
    expr = foreach(fl = files, .combine = merge_all) %dopar% {
        samplemap_table(fl, soft, GSE)
    }
    warnings()
    if (GSE %in% names(special_resolver)) {
        also = special_resolver[[GSE]](soft)
        expr = merge_all(expr, also)
    }
    samps = colnames(expr)
    dups = samps[duplicated(samps)]
    no_int = samps[!is.int(expr)]
    unusful = intersect(dups, no_int)
    ind = colnames(expr) %in% unusful
    expr = expr[, !ind]
    colnames(expr) = make.unique(colnames(expr))
    return(expr)
}

is.int = function(x, MARGIN = 2, na.rm = FALSE) {
    return(apply(x%%1==0, MARGIN, all, na.rm = na.rm))
}

get_gsm_title = function(soft) {
    return(sapply(soft@gsms, function(x) x@header$title))
}

get_gsm_field = function(soft, f_name, line_number = T) {
    return(sapply(soft@gsms, function(x)
                  x@header[[f_name]][line_number]))
}

samplemap_table = function(fl, soft, GSE) {
    if (basename(fl) %in% exclude) return(NULL)
    txt = readLines(fl, 10)
    comments = grepl("^[#%]", txt)
    if (all(comments)) {
        warning(paste("No counts extracted for", fl))
        return(NULL)
    }
    skip = which(!comments)[1]
    header = txt[skip]
    samps = strsplit(header, "\t")[[1]]
    if (length(samps) == 1) {
        samps = strsplit(header, ",")[[1]]
        if (length(samps) == 1) {
            samps = strsplit(header, ";")[[1]]
            expr = read.csv2(fl, header = FALSE, skip = skip,
                             stringsAsFactors = FALSE)
        } else {
            expr = read.csv(fl, header = FALSE, skip = skip,
                            stringsAsFactors = FALSE)
        }
    } else {
        expr = read.delim(fl, header = FALSE, skip = skip,
                          stringsAsFactors = FALSE)
    }
    rownames(expr) = make.unique(expr[, 1])
    expr = expr[, -1, drop = FALSE]
    lastCol = ncol(expr)
    if (all(is.na(expr[lastCol]))) expr = expr[-lastCol]
    if (GSE %in% names(samplename_resolver)) {
        samps = samplename_resolver[[GSE]](fl, soft)
        if (is.null(samps)) return(NULL)
        samps = rep(samps, length.out = ncol(expr)) # ensure length
        ind = !is.na(samps)
        expr = expr[, ind]
        colnames(expr) = make.unique(samps[ind])
        return(expr)
    }
    if (grepl("^GSM[0-9]", basename(fl))) {
        # Sinle cell experiment
        number = apply(expr, 2, is.numeric)
        if (any(!number) || any(expr%%1!=0)) {
            warning(paste("Are thouse single cell counts:", fl))
            return(NULL)
        }
        samps = gsub("(GSM[0-9]*).*", "\\1", basename(fl))
        samps = rep(samps, length.out = ncol(expr))
        colnames(expr) = make.unique(samps)
        return(expr)
    }
    samps = samps[samps != ""]
    if (length(samps) == ncol(expr)+1) {
        samps = samps[-1]
    } else if (length(samps) == ncol(expr)+2) {
        samps = samps[samps != ""]
    }
    if (length(samps) != ncol(expr)) {
        warning(paste("Could not match header of file", fl))
    }
    return(generic_map(expr, soft, samps, name=fl))
}

generic_map = function(expr, soft, samps, name=NULL) {
    titles = sapply(soft@gsms, function(x) x@header$title)
    if (length(samps) > length(titles)) {
        tol = length(samps) - length(titles)
    } else {
        tol = length(intersect(samps, not_a_sample))
    }
    map = match(samps, titles)
    if (sum(is.na(map)) > tol) {
        nospaceIn = toupper(gsub("[ -._]", "", samps))
        nospaceTitle = toupper(gsub("[ -._]", "", titles))
        map = pmatch(nospaceIn, nospaceTitle, dup = TRUE)
    }
    for (i in 1:5) {
        if (sum(is.na(map)) > tol) {
            descs = get_gsm_field(soft, "description", i)
            descs = gsub("^[^:]*: ", "", descs)
            map = pmatch(samps, descs, dup = TRUE)
        }
        if (sum(is.na(map)) > tol) {
            tryCatch({
                chr = get_gsm_field(soft, "characteristics_ch1", i)
                chr = gsub("^[^:]*: ", "", chr)
                map = pmatch(samps, chr)
            })
        }
    }
    if (sum(is.na(map)) > tol) {
        warning(paste("Sample in table", name, "could not be matched."))
        return(NULL)
    }
    colnames(expr) = names(soft@gsms)[map]
    expr = expr[, !is.na(map)]
    return(expr)
}

custom_extract = function(GSE){
    return(custom_resolvers[[GSE]]())
}

mmseq_extract = function(files) {
    if (length(files) == 0) return(NULL)
    expr = foreach(fl = files, .combine = merge_all) %dopar% {
        res = read.table(fl, header = TRUE)
        rownames(res) = res$feature_id
        return(res[, "unique_hits", drop = FALSE])
    }
    gsms = gsub("_.*", "", basename(files))
    colnames(expr) = gsms
    return(expr)
}

kallisto_extract = function(files) {
    if (length(files) == 0) return(NULL)
    expr = foreach(fl = files, .combine = merge_all) %dopar% {
        res = read.table(fl, header = TRUE)
        if ("gene_id" %in% colnames(res)) {
            rownames(res) = res$gene_id
            res = res[, "expected_count", drop = FALSE]
        } else {
            rownames(res) = res$target_id
            res = res[, "est_counts", drop = FALSE]
        }
        return(floor(res))
    }
    gsms = gsub("_.*", "", basename(files))
    colnames(expr) = gsms
    return(expr)
}

excel_extract = function(files) {
    if (length(files) == 0) return(NULL)
    GSE = basename(dirname(files[1]))
    soft = get_soft(GSE)
    expr = foreach(fl = files, .combine = merge_all) %dopar% {
        res = read.xlsx(fl)
        rownames(res) = make.unique(res[, 1])
        res = res[, -1, drop = FALSE]
        samps = colnames(res)
        return(generic_map(res, soft, samps, name=fl))
    }
    return(expr)
}
