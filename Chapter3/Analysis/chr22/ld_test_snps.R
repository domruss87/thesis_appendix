args = commandArgs(trailingOnly = TRUE)
args_file = args[1]
bim_loc = args[2]
if (length(args) == 3){
    prefix = args[3]
} else {
    prefix = ""
}

#args = "off_jr2_aaBb_or1.2_maf0.1_it1"

library(tidyverse)
library(data.table)

sd = as.integer(substr(gsub("[^0-9]","",args_file),1,6))
set.seed(sd)

blk = fread(paste0("ld_test/",args_file,".blocks"), sep = "\t", header = FALSE, col.names = "snps")

blk$snps = gsub("^\\* ","",blk$snps)

snplist = c()
for (i in 1:nrow(blk)){
	addsnp = sample(unlist(strsplit(as.character(blk[i,]), split = " ")), 1)
    if(grepl("rs_int1", as.character(blk[i,])) & !grepl("rs_int1", addsnp)){
        addsnp = paste0("rs_int1p_",addsnp)
    } else if(grepl("rs_int2", as.character(blk[i,])) & !grepl("rs_int2", addsnp)){
        addsnp = paste0("rs_int2p_",addsnp)
    }
	snplist = c(snplist, addsnp)
}

blk_snps = unlist(strsplit(as.character(blk$snps), split = " "))

bim = fread(bim_loc, col.names = c("chr","id","cm","bp","a1","a2"), header = FALSE)

outstanding = bim$id[!bim$id %in% blk_snps]

analyze = data.frame(id = c(snplist,outstanding))

write_tsv(analyze, paste0("ld_test/",prefix,args_file,"_snps.txt"), col_names = FALSE)
