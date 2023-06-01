args = commandArgs(trailingOnly = TRUE)
library(data.table)
library(tidyverse)

dat = fread(args, sep = " ")
print("data in")


snps = colnames(dat)[2:length(colnames(dat))]

dat = dat %>% select(all_of(snps),PHENOTYPE)

dat$PHENOTYPE[dat$PHENOTYPE == 1] = 0
dat$PHENOTYPE[dat$PHENOTYPE == 2] = 1

write_tsv(dat, args)
