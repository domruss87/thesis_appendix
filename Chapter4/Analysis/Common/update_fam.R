library(tidyverse)
library(data.table)

fam = fread("plink/filtered_imputed.fam", col.names = c("fid","iid","pid","mid","sex","p"))
ped = fread("DDD10K_family_relationships.txt") %>% filter(individual_id %in% fam$iid & dad_id != 0 )

for (i in fam$iid){
	if(i %in% ped$individual_id){
		pat = ped$dad_id[ped$individual_id == i]
		mat = ped$mum_id[ped$individual_id == i]
		sex = ped$sex[ped$individual_id == i]
		if(sex == "M"){
			sex = 1
		} else {
			sex = 2
		}
		
		fam$pid[fam$iid == i] = pat
		fam$mid[fam$iid == i] = mat
		fam$sex[fam$iid == i] = sex
		fam$sex[fam$iid == pat] = 1
		fam$sex[fam$iid == mat] = 2
		fam$p[fam$iid == i] = 2
		fam$p[fam$iid == pat] = 1
		fam$p[fam$iid == mat] = 1
		
		fam$fid[fam$iid %in% c(i,pat,mat)] = paste0("f",i)
	}
}

fam$fid[fam$iid %in% c("DDDP100248","DDDP110977","DDDP107066","DDDP107067")] = paste0("fDDDP100248","_","fDDDP110977")
fam$fid[fam$iid %in% c("DDDP101540","DDDP101541","DDDP101539","DDDP101538")] = paste0("fDDDP101540","_","fDDDP101541")


write_tsv(fam, "plink/filtered_imputed.fam", col_names = F)
