library(data.table)

file="DREAM_r3_coarse.csv"
dat = read.csv(file)

sdat = subset(dat, metric=="spearman")
dt = data.table(sdat)
mms = dt[, list(mean=mean(metric_value), median=median(metric_value)),
         by=c("repo_name", "submitter", "is_latest")]
mdt1 = data.frame(mms[order(mean, decreasing = TRUE)])
mdt2 = data.frame(mms[order(median, decreasing = TRUE)])

print(mdt1)
