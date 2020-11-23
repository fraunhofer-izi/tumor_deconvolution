
RNASeq_file = "../data/extraction_stats.RDS"
array_file = "../data/extraction_stats_array.RDS"
out_file_base = "../data/annotation"

sdf = readRDS(RNASeq_file)
adf = readRDS(array_file)

sdf$experiment = "rna-seq"
adf$experiment = "expression-array"

cols = union(colnames(sdf), colnames(adf))

for (co in cols) {
    if (!co %in% colnames(sdf)) {
        sdf[, co] = NA
    }
    if (!co %in% colnames(adf)) {
        adf[, co] = NA
    }
}

map1 = match(cols, colnames(sdf))
map2 = match(cols, colnames(adf))

ndf = rbind(sdf[, map1], adf[, map2])


saveRDS(ndf, paste0(out_file_base, ".RDS"))
write.table(ndf, paste0(out_file_base, ".tsv"), quote = FALSE,
            fileEncoding = "utf-8", sep="\t")
