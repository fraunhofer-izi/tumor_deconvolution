library(foreach)
library(doMC)
registerDoMC(32)

outPath = "../data/"
outFile = file.path(outPath, "HUGO_mpa.RDS")
HUGOfile = file.path(outPath, "HUGO.tsv")
HUGO_df = read.delim(HUGOfile, colClasses = "character")
hugo_col = HUGO_df$Approved.symbol
HUGO_map = foreach(i = 1:ncol(HUGO_df), .combine = rbind) %do% {
    message(paste("column", colnames(HUGO_df)[i]))
    index = HUGO_df[, i] != ""
    ref_name = colnames(HUGO_df)[i]
    map = data.frame(HUGO = hugo_col[index],
                     ref = toupper(HUGO_df[index, i]),
                     source = ref_name,
                     stringsAsFactors = FALSE)
    broken = grepl(", ", map$ref)
    if (any(broken)) {
        fixed = strsplit(map$ref[broken], ", ")
        broken_ind = which(broken)
        fixed = foreach(j = 1:length(broken_ind), .multicombine = TRUE,
                        .combine = rbind) %dopar% {
            data.frame(HUGO = map$HUGO[broken_ind[j]],
                       ref = fixed[[j]],
                       source = ref_name)
        }
        map = rbind(subset(map, !broken), fixed)
    }
    return(map)
}
HUGO_map = unique(HUGO_map)
map_count = table(HUGO_map$ref)
HUGO_map$matches = map_count[HUGO_map$ref]
no_HUGO = c(HUGO = "none", ref = NA, source = "none", matchaes = 1)
HUGO_map = rbind(no_HUGO, HUGO_map)
HUGO_map$HUGO = as.factor(HUGO_map$HUGO)

saveRDS(HUGO_map, outFile)
