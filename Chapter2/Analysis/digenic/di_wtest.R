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

res = wtest(geno, pheno, w.order = 2, hf1 = "default.hf1",hf2 = "default.hf2", which.marker = NULL, output.pval = NULL,sort = TRUE, input.pval = 1, input.poolsize = NULL)

write_tsv(res$results, paste0(output,"_wtest_2nd_ord.tsv"))

saveRDS(res, paste0(output,"_wtest_2nd_ord.rds"))
