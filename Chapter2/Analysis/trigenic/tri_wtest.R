args = commandArgs(trailingOnly = TRUE)
input = args[1]
output = args[2]

library(wtest)
library(tidyverse)
library(data.table)

fin = fread(paste0(input,".raw"), sep = " ", header = T)

geno = fin %>% select(7:ncol(fin))
pheno = fin$PHENOTYPE
pheno[pheno == 1] = 0
pheno[pheno == 2] = 1

step1start = timestamp()

res.high = wtest.high(geno, pheno, w.order = 3, hf1 = "default.hf1", hf.high.order = "default.high", which.marker = NULL, output.pval = NULL, sort = TRUE, input.pval = NULL, input.poolsize = NULL)

step2done <- timestamp()
iRFtimes <- c(step1start, step2done)
rnames <- c("Start", "Complete")
mylog <- data.frame(iRFtimes, row.names = rnames)
write_tsv(mylog, paste0(output,"_wtest_timings.tsv"))
write_tsv(res.high$results, paste0(output,"_wtest_3rd_ord.tsv"))

saveRDS(res.high, paste0(output,"_wtest_3rd_ord.rds"))
