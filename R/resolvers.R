#! /usr/bin/env Rscript
suppressPackageStartupMessages(library(foreach))

source("R/utilities.R")

custom_resolvers = list(
    "GSE116555" = function(GSE = "GSE116555") {
        file_name = "GSE116555_combine_alignment.Rda"
        file_path = file.path(studies, GSE, file_name)
        data_env = attemp_load(file_path)
        expr = data.frame(dat@raw.data)

        # sample names
        samps = colnames(expr)
        d8 = grepl("^D8", samps)
        samps[d8] = "GSM3242729"
        d12 = grepl("^D12", samps)
        samps[d12] = "GSM3242730"
        colnames(expr) = make.unique(samps)
        expr = expr[, d8 | d12]
        return(expr)
    },
    "GSE114519" = function(GSE = "GSE114519") {
        file_name = "GSE114519_RawCounts.txt"
        file_path = file.path(studies, GSE, file_name)
        expr = read.table(file_path, header = TRUE, row.names = 1)
        expr = expr[, -1]

        # sample names
        samps = colnames(expr)
        soft = get_soft(GSE)
        titles = get_gsm_title(soft)
        titles = make.names(titles)
        map = match(samps, titles)
        colnames(expr) = names(soft@gsms)[map]
        return(expr)
    },
    "GSE92506" = function(GSE = "GSE92506") {
        file_name = "GSE92506_RawCounts.txt"
        file_path = file.path(studies, GSE, file_name)
        expr = read.table(file_path, header = TRUE, row.names = 1)
        expr = expr[, -1]

        # sample names
        samps = colnames(expr)
        soft = get_soft(GSE)
        titles = get_gsm_title(soft)
        map = match(samps, titles)
        colnames(expr) = names(soft@gsms)[map]
        return(expr)
    },
    "GSE107011" = function(GSE = "GSE107011") {
        file_name = "GSE107011_Processed_data_TPM.txt"
        file_path = file.path(studies, GSE, file_name)
        expr = read.table(file_path)

        # sample names
        samps = colnames(expr)
        samps = gsub("\\.", "p", samps)
        samps = gsub("pp1", "m", samps)
        soft = get_soft(GSE)
        descs = sapply(soft@gsms, function(gsm) {
                           tit = gsm@header$title
                           res = sub("_rep[0-9]*$", "", tit)
                           res = gsub("\\+$", "p", res) 
                           res = gsub("-$", "m", res) 
                           return(res)
        })
        descs = make.names(descs)
        descs = gsub("\\.", "p", descs)
        map = match(samps, descs)
        colnames(expr) = names(soft@gsms)[map]
        return(expr)
    },
    "GSE64655" = function(GSE = "GSE64655") {
        file_name = paste0("GSE64655_Normalized_transcript_expression_in",
                           "_human_immune_cells.txt")
        file_path = file.path(studies, GSE, file_name)
        expr = read.delim(file_path, skip = 4, header = TRUE)

        # feature names
        rownames(expr) = expr$Gene.ID

        # sample names
        samps = colnames(expr)
        soft = get_soft(GSE)
        donor = get_gsm_field(soft, "characteristics_ch1", 4)
        donor = gsub(".*: ", "", donor)
        ctype = sct = get_gsm_field(soft, "characteristics_ch1", 2)
        sct[grepl("PBMC$", ctype)] = "PBMC"
        sct[grepl("myeloid DC$", ctype)] = "mDC"
        sct[grepl("monocytes$", ctype)] = "Mono"
        sct[grepl("neutrophils$", ctype)] = "Neut"
        sct[grepl("B cells$", ctype)] = "B"
        sct[grepl("T cells$", ctype)] = "T"
        sct[grepl("NK cells$", ctype)] = "NK"
        day = get_gsm_field(soft, "characteristics_ch1", 1)
        day = gsub(".*: ", "", day)
        day = gsub(" d", "", day)
        day = paste0("d", day)
        lables = paste(donor, sct, day, sep = "_")
        map = match(samps, lables)
        colnames(expr) = names(soft@gsms)[map]
        expr = expr[, !is.na(map)]
        return(expr)
    },
    "GSE49159" = function(GSE = "GSE49159") {
        file_name = "GSE49159_RPKM_values.txt"
        file_path = file.path(studies, GSE, file_name)
        expr = read.delim(file_path, header = TRUE)

        # feature names
        rownames(expr) = expr$geneName

        # sample names
        samps = colnames(expr)
        samps = gsub("Count_", "GLD-", samps)
        soft = get_soft(GSE)
        titles = get_gsm_title(soft)
        map = match(samps, titles)
        colnames(expr) = names(soft@gsms)[map]
        expr = expr[, !is.na(map)]
        return(expr)
    },
    "GSE107019" = function(GSE = "GSE107019") {
        # could only find the microarray data
        return(NULL)
    },
    "GSE118165" = function(GSE = "GSE118165") {
        file_name = "GSE118165_RNA_gene_abundance.txt"
        file_path = file.path(studies, GSE, file_name)
        expr = read.table(file_path)

        # sample names
        samps = colnames(expr)
        soft = get_soft(GSE)
        descs = sapply(soft@gsms, function(gsm) {
                           return(gsm@header$title)
        })
        descs = make.names(descs)
        map = match(samps, descs)
        colnames(expr) = names(soft@gsms)[map]
        return(expr)
    },
    "GSE114135" = function(GSE = "GSE114135") {
        # no seq data found
        return(NULL)
    }
)

samplename_resolver = list(
    "GSE103568" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), "\t")[[1]]
        samps = paste0("UHRR_", samps)
        titles = get_gsm_title(soft)
        map = match(samps, titles)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE100576" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), "\t")[[1]]
        titles = get_gsm_title(soft)
        t1 = t2 = ifelse(grepl("Mock", titles), "mock", "WT")
        ind = grepl("^Cytoplasmic", titles)
        t2[ind] = "Cy"
        ind = grepl("^Total", titles)
        t2[ind] = "To"
        ind = grepl("^Chromatin-associated", titles)
        t2[ind] = "Cr"
        ind = grepl("^Nucleoplasmic", titles)
        t2[ind] = "Nu"
        t3 = ifelse(grepl("1$", titles), "1", "2")
        short_title = paste(t1, t2, t3, sep = "_")
        map = match(samps, short_title)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE100576" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), "\t")[[1]]
        titles = get_gsm_title(soft)
        t1 = ifelse(grepl("^Healthy", titles), "N", "D")
        t2 = gsub("^[^_]*_([^_]*)_.*$", "\\1", titles)
        t3 = gsub(".*([0-9])$", "\\1", titles)
        short_title = paste(t1, t2, t3, sep = "-")
        map = match(samps, short_title)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE57353" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), "\t")[[1]][-1]
        samps = toupper(samps)
        samps = gsub("_", "", samps)
        titles = get_gsm_title(soft)
        titles = toupper(titles)
        titles = gsub(" ", "", titles)
        titles = gsub("CONTROL1", "CONTROL", titles)
        map = pmatch(samps, titles)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE74324" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), "\t")[[1]]
        titles = get_gsm_title(soft)
        res = rep(NA, length.out = length(samps))
        res[1:length(titles)] = names(soft@gsms)
        return(res)
    },
    "GSE77443" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), "\t")[[1]][-1]
        samps = toupper(samps)
        samps = gsub("_", "", samps)
        samps = gsub("RPKM", "", samps)
        titles = get_gsm_title(soft)
        titles = gsub("Replicate", "", titles)
        titles = toupper(titles)
        titles = gsub(" ", "", titles)
        map = pmatch(samps, titles)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE111880" = function(fl, soft) {
        res = names(soft@gsms)
        return(res)
    },
    "GSE90569" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), "\t")[[1]][-1]
        samps = gsub("h_", "_", samps)
        titles = get_gsm_title(soft)
        map = pmatch(samps, titles)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE106542" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), "\t")[[1]][-1]
        samps = gsub("^BRNA_", "", samps)
        titles = get_gsm_title(soft)
        titles = gsub("Bulk_RNA-seq_", "", titles)
        map = pmatch(samps, titles)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE109843" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), ",")[[1]][-1]
        samps = gsub("_.*", "", samps)
        if (basename(fl) == "GSE109843_project1.csv") i = 3
        if (basename(fl) == "GSE109843_project2.csv") i = 2
        descs = get_gsm_field(soft, "description", i)
        map = pmatch(samps, descs)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE41166" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), "\t")[[1]][-1]
        descs = get_gsm_field(soft, "description", 5)
        map = pmatch(samps, descs)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE93672" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), "\t")[[1]]
        samps = paste("RNA-seq", samps)
        titles = get_gsm_title(soft)
        map = pmatch(samps, titles)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE109448" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), ",")[[1]][-1]
        samps = paste("RNA-seq", samps)
        titles = get_gsm_title(soft)
        map = pmatch(samps, titles)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE102045" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), ",")[[1]][-1]
        samps = gsub("^JM[0-9]*-", "", samps)
        samps = gsub("-[0-9]*$", "", samps)
        titles = get_gsm_title(soft)
        map = pmatch(samps, titles)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE102250" = function(fl, soft) {
        if (basename(fl) == "GSE102250_counts.txt") {
            # srange sample names while all are in other file
            return(NULL)
        }
        samps = strsplit(readLines(fl, 1), ",")[[1]][-1]
        titles = get_gsm_title(soft)
        titles = gsub("control", "ctrl", titles)
        titles = gsub(" ", "_", titles)
        map = pmatch(samps, titles)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE70330" = function(fl, soft) {
        titles = get_gsm_title(soft)
        res = names(soft@gsms)
        return(res)
    },
    "GSE77108" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), "\t")[[1]][-1]
        samps = gsub("-", "_", samps)
        samps = gsub("_([0-9])$", "_R\\1", samps)
        titles = get_gsm_title(soft)
        titles = gsub("Healthy", "N", titles)
        titles = gsub("Diabetic", "D", titles)
        titles = gsub("_siRNA", "", titles)
        map = pmatch(samps, titles)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE100562" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), "\t")[[1]]
        samps = gsub("[^_]*_", "", samps)
        libs = get_gsm_field(soft, "characteristics_ch1", 4)
        libs = gsub("ID: ", "", libs)
        map = match(samps, libs)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE85530" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), "\t")[[1]][-1]
        samps = gsub("[^_]*_", "", samps)
        libs = get_gsm_field(soft, "characteristics_ch1", 1)
        libs = gsub("samplelabel: ", "", libs)
        map = match(samps, libs)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE111907" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), "\t")[[1]][-1]
        samps = gsub("_.*$", "", samps)
        libs = get_gsm_field(soft, "characteristics_ch1", 1)
        libs = gsub("patient id: ", "", libs)
        map = match(samps, libs)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE96975" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), "\t")[[1]]
        samps = samps[samps != ""]
        samps = toupper(samps)
        samps = gsub("_", "", samps)
        titles = get_gsm_title(soft)
        titles = gsub("_1$", "", titles)
        titles = gsub("_", "", titles)
        titles = toupper(titles)
        map = pmatch(samps, titles)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE87186" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), "\t")[[1]][-1]
        days = get_gsm_field(soft, "characteristics_ch1", 3)
        days = gsub("day: ", "", days)
        ids = get_gsm_field(soft, "characteristics_ch1", 4)
        ids = gsub("subjectid: ", "", ids)
        lables = paste0(ids, "_D", days)
        map = match(samps, lables)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE90600" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), "\t")[[1]][-1]
        descs = get_gsm_field(soft, "description")
        map = match(samps, descs)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE85112" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), "\t")[[1]]
        n = length(samps) - 4
        samps = gsub("_.*", "", basename(fl))
        samps = rep(samps, length.out = n)
        res = c(NA, NA, NA, samps)
        return(res)
    },
    "GSE90711" = function(fl, soft) {
        res = names(soft@gsms)
        return(res)
    },
    "GSE109449" = function(fl, soft) {
        if (basename(fl) != "singlecell_rnaseq_gene_counts.tsv") {
            return(NULL)
        }
        samps = strsplit(readLines(fl, 1), "\t")[[1]][-1]
        titles = get_gsm_title(soft)
        titles = gsub(".*\\((.*)\\).*", "\\1", titles)
        map = match(samps, titles)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE65515" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), "\t")[[1]][-1]
        samps = gsub("_[^_]*$", "", samps)
        titles = get_gsm_title(soft)
        ind = grepl("RNA-seq", titles) & ! grepl("small", titles)
        titles = titles[ind]
        titles = gsub(" .*", "", titles)
        titles = gsub("-[0-9]*", "", titles)
        map = match(samps, titles)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE60996" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), "\t")[[1]][-1]
        titles = get_gsm_title(soft)
        titles = gsub("_iPSC", "", titles)
        map = match(samps, titles)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE100745" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), ",")[[1]][-1]
        samps = gsub("Naive_AKJ_D3630", "Naive_AKJ", samps)
        samps = gsub("Naive_EVC422_D4745", "Naive_EVC422", samps)
        titles = get_gsm_title(soft)
        map = match(samps, titles)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE115240" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), "\t")[[1]][-1]
        samps = gsub("_C", "_O", samps) # typo ?
        titles = get_gsm_title(soft)
        map = match(samps, titles)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE42212" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), "\t")[[1]][-1]
        samps = gsub("\\(RPKM\\)", "", samps)
        titles = get_gsm_title(soft)
        titles = gsub("_RNAseq", "", titles)
        map = match(samps, titles)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE119088" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), "\t")[[1]][-1]
        titles = get_gsm_title(soft)
        map = match(samps, titles)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE50781" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), "\t")[[1]][-1]
        titles = get_gsm_title(soft)
        map = match(samps, sort(samps))
        res = sort(names(titles))[map]
        return(res)
    },
    "GSE78922" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), "\t")[[1]][-1]
        samps = gsub("-.*", "", samps)
        descs = get_gsm_field(soft, "description")
        map = match(samps, descs)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE54384" = function(fl, soft) {
        titles = get_gsm_title(soft)
        others = rep(NA, 6)
        res = c(names(titles), others)
        return(res)
    },
    "GSE98758" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), "\t")[[1]][-1]
        titles = get_gsm_title(soft)
        map = match(samps, titles)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE106540" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), "\t")[[1]][-1]
        samps = gsub("SCRNA_", "", samps)
        titles = get_gsm_title(soft)
        titles = gsub("Single-cell_RNA-seq_", "", titles)
        map = match(samps, titles)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE75643" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), "\t")[[1]][-1]
        samps = make.names(samps)
        descs = get_gsm_field(soft, "description")
        map = match(samps, descs)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE92532" = function(fl, soft) {
        samps = strsplit(readLines(fl, 1), "\t")[[1]][-1]
        descs = get_gsm_field(soft, "description", 3)
        if (grepl("batch1", fl)) {
            samps = gsub("-", "_", samps)
            samps = paste0("S", samps)
            descs = gsub("\\.", "_", descs)
        }
        map = match(samps, descs)
        res = names(soft@gsms)[map]
        return(res)
    },
    "GSE80361" = function(fl, soft) {
        samp = gsub(".*/(GSM[0-9]*)_[^/]*", "\\1", fl)
        res = c(samp, NA)
        return(res)
    },
    "GSE101341" = function(fl, soft) {
        if (grepl("CD34_umis.txt$", fl)) {
            return("GSM2701545")
        }
        design_path = file.path(dirname(fl), "GSE101341_experimental_design.txt")
        desig = read.delim(design_path)
        samps = strsplit(readLines(fl, 1), "\t")[[1]][-1]
        dmap = match(samps, desig$Well_ID)
        tit = paste("CD14+", desig$DPI[dmap], "dpi")
        titles = get_gsm_title(soft)
        map = match(tit, titles)
        return(names(titles)[map])
    }
)

special_resolver = list(
    "GSE101341" = function(soft) {
        lpath = file.path(studies, "GSE101341")
        fl = file.path(lpath, "GSE101341_CD14_umis.txt")
        studie = file.path(lpath, "GSE101341_experimental_design.txt")
        studie = read.delim(studie)
        expr = read.delim(fl, row.names = 1, header = FALSE, skip = 1)
        samps = strsplit(readLines(fl, 1), "\t")[[1]][-1]
        map = match(samps, studie$Well_ID)
        dpi = studie$DPI[map]
        samp_t = paste("CD14+", dpi, "dpi")
        titles = get_gsm_title(soft)
        map = match(samp_t, titles)
        samps = names(soft@gsms)[map]
        colnames(expr) = samps
        return(expr)
    }
)

filter_additionals = list(
    "GSE72056" = function(expr, cell_type) {
        cell_type = make.names(cell_type)
        map = c(T.cells = 1,
                B.cells = 2,
                macrophages = 3,
                monocytic.lineage = 3,
                endothelial.cells = 4,
                fibroblasts = 5,
                NK.cells = 6)
        if (! cell_type %in% names(map)) return(NULL)
        indrow = "non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)"
        ind = expr[indrow, ] == map[cell_type]
        return(expr[4:nrow(expr), ind, drop = FALSE])
    }
)

# files and studies with no usful data
exclude = c("GSE104174_TPM_SSc_study.txt",
            "GSE62056_EZH2_Karpas422_RNA-seq_summary.txt",
            "GSE98758_humanPromoterChIPcounts.txt",
            "GSE56788", "GSE41687", "GSE100250", "GSE112525",
            "GSE85573")
