args = commandArgs(trailingOnly = TRUE)

library(tidyverse)
library(data.table)
library(snpStats)

#sd = as.integer(substr(gsub("[^0-9]","",args_file),1,6))
#set.seed(sd)
set.seed(123)
print(1)
#load genetic data
gen = fread("chr22.raw") 
gen = gen %>% 
  dplyr::select(-FID,-IID,-PAT,-MAT,-SEX,-PHENOTYPE) 

colnames(gen) = gsub("_[A-Z]*$","",colnames(gen))
colnames(gen) = gsub("Affx-","rs999x",colnames(gen))

gen$pheno = sample(c(rep(0,50000),rep(1,50000)), 100000, replace = FALSE)

#positions and alleles
bim = fread("chr22.bim", col.names = c("chr","id","cm","bp","a1","a2"))

bim$id = gsub("Affx-","rs999x",bim$id)

#ld data
ld_lkup = fread("chr22.ld.gz")
print(2)
#snp summaries
snp_info_func = function(snp){
  maf = sum(snp, na.rm = T)/(2*length(snp))
  num2 = length(snp[snp == 2 & !is.na(snp)])
  num1 = length(snp[snp == 1 & !is.na(snp)])
  num0 = length(snp[snp == 0 & !is.na(snp)])
  lout = c(maf,num2,num1,num0)
  return(lout)
}

maf_func = function(snp){
  maf = sum(snp, na.rm = T)/(2*length(snp))
  return(maf)
}

snp_info = matrix(ncol = 10, nrow = 0)

mafs = lapply(gen, maf_func)
case_snps = lapply(gen[gen$pheno == 1], snp_info_func)
case_snps = do.call("rbind",case_snps)

cont_snps = lapply(gen[gen$pheno == 0], snp_info_func)
cont_snps = do.call("rbind",cont_snps)

snp_info = as.data.frame(cbind(case_snps, cont_snps))
colnames(snp_info) = c("case_maf","case_2","case_1","case_0","cont_maf","cont_2","cont_1","cont_0")
snp_info$maf = mafs
print(3)
#example
dl_loc = args
mod_dat = fread(dl_loc) %>% select("M0P0","M0P1","Class")
print(head(mod_dat))

mod_info1_case = lapply(mod_dat[mod_dat$Class==1,], snp_info_func)[[1]]
mod_info1_cont = lapply(mod_dat[mod_dat$Class==0,], snp_info_func)[[1]]
mod_info2_case = lapply(mod_dat[mod_dat$Class==1,], snp_info_func)[[2]]
mod_info2_cont = lapply(mod_dat[mod_dat$Class==0,], snp_info_func)[[2]]

print(mod_info1_case)
print(mod_info2_case)
print(mod_info2_case)
print(mod_info2_cont)

ln_cnt = 0
bnd = 0.01
while(ln_cnt < 50){ 
    print(bnd)
    int_options = snp_info %>% filter(case_2 > mod_info1_case[2] - bnd * mod_info1_case[2] & 
               case_2 < mod_info1_case[2] + bnd * mod_info1_case[2] &
               cont_2 > mod_info1_cont[2] - bnd * mod_info1_cont[2] & 
               cont_2 < mod_info1_cont[2] + bnd * mod_info1_cont[2])
    ln_cnt = nrow(int_options)
    bnd = bnd + 0.01
}

mod22_case = mod_dat %>% filter(Class == 1) %>% dplyr::select(-Class) %>% filter_all(all_vars(. == 2)) %>% nrow
mod22_cont = mod_dat %>% filter(Class == 0) %>% dplyr::select(-Class) %>% filter_all(all_vars(. == 2)) %>% nrow
mod22 = mod22_case + mod22_cont
print(4)
df = expand.grid(a = rownames(int_options), b = rownames(int_options))
print(df)
df = df %>%
  filter(a != b) %>%
  rowwise() %>% 
  mutate(a = as.character(a),
         b = as.character(b),
         int = paste(sort(c(a, b)), collapse = "_")) %>%  
  dplyr::select(-a,-b) %>%
  unique %>%
  separate(int, into = c("a","b"), sep = "_")

print(df)
print(nrow(df))

ld_here = ld_lkup %>%
  filter(SNP_A %in% rownames(int_options) & SNP_B %in% rownames(int_options))

diffs = c()
print(5)
for (i in 1:nrow(df)){
  snp1 = df[[i,1]]
  snp2 = df[[i,2]]
  print(5.0)
  joint_ld = ld_here %>% filter((SNP_A == snp1 & SNP_B == snp2)|(SNP_A == snp2 & SNP_B == snp1)) %>% nrow
  print(5.1)
  if(joint_ld == 0){
    num22_case = gen %>%
      dplyr::select(pheno, all_of(snp1), all_of(snp2)) %>%
      filter(pheno == 1) %>%
      filter_at(vars(starts_with("rs")), all_vars(. == 2)) %>%
      nrow
  print(5.2)
    num22_cont = gen %>%
      dplyr::select(pheno, all_of(snp1), all_of(snp2)) %>%
      filter(pheno == 0) %>%
      filter_at(vars(starts_with("rs")), all_vars(. == 2)) %>%
    nrow
   print(5.3)   
    d = abs(mod22 - (num22_case + num22_cont))
    diffs = c(diffs, d)
    int_name = paste0(snp1,"_",snp2)
    names(diffs)[length(diffs)] = int_name
  print(5.4)
  }
}
print(5.5)

best_pair = names(diffs[diffs == min(diffs)][1])
s1 = unlist(strsplit(best_pair, split = "_"))[1]
s2 = unlist(strsplit(best_pair, split = "_"))[2]

#now build pair next to these
gen = gen %>%
  add_column(rs_new_var1 = NA, .after = s1) %>%
  add_column(rs_new_var2 = NA, .after = s2) 
gen$rs_new_var1 = as.integer(gen$rs_new_var1)
gen$rs_new_var2 = as.integer(gen$rs_new_var2)

colnames(gen)[match(s1, colnames(gen))] = "s1"
colnames(gen)[match(s2, colnames(gen))] = "s2"

#table(paste(int_gen$rs_new_var1,int_gen$rs_new_var2,int_gen$pheno, sep = "_"))
#table(paste(mod_dat$M0P0,mod_dat$M0P1, mod_dat$Class, sep = "_"))
#gen = gen2
#looping
gt_sorted = names(sort(table(paste(mod_dat$M0P0,mod_dat$M0P1,mod_dat$Class, sep = "_"))))

gen$id = paste0("id",1:nrow(gen))
int_gen = gen %>% dplyr::select(id,s1,s2,rs_new_var1,rs_new_var2,pheno)
int_gen$rs_new_var1[is.na(int_gen$s1)] = 5
int_gen$rs_new_var2[is.na(int_gen$s2)] = 5
for (gt in gt_sorted){
  print(gt)
  g1 = strsplit(gt, split = "_")[[1]][1]
  g2 = strsplit(gt, split = "_")[[1]][2]
  p1 = strsplit(gt, split = "_")[[1]][3]
  
  gtr = int_gen %>% filter(s1 == g1 & s2 == g2 & pheno == p1 & 
                             (is.na(rs_new_var1) | rs_new_var1 == 5) & (is.na(rs_new_var2) | rs_new_var2 == 5) ) %>% nrow
  gtm = mod_dat %>% filter(M0P0 == g1, M0P1 == g2, Class == p1) %>% nrow

 print(table(paste(int_gen$rs_new_var1,int_gen$rs_new_var2,int_gen$pheno, sep = "_")))
print(table(paste(mod_dat$M0P0,mod_dat$M0P1, mod_dat$Class, sep = "_"))) 
print(6)
  if(gtr <= gtm){
    int_gen$rs_new_var1[int_gen$s1 == g1 & 
                        int_gen$s2 == g2 & 
                        int_gen$pheno == p1 &
                        (is.na(int_gen$rs_new_var1) | int_gen$rs_new_var1 == 5) & 
                        (is.na(int_gen$rs_new_var2) | int_gen$rs_new_var2 == 5)] = g1
                    
    int_gen$rs_new_var2[int_gen$s1 == g1 & 
                        int_gen$s2 == g2 &
                        int_gen$pheno == p1 &
                        int_gen$rs_new_var1 == g1 & 
                        (is.na(int_gen$rs_new_var2) | int_gen$rs_new_var2 == 5)] = g2
print(6.5)
                    
#    gtn = int_gen %>% filter(rs_new_var1 == g1 & rs_new_var2 == g2 & pheno == p1) %>% nrow
#    gtm = mod_dat %>% filter(M0P0 == g1, M0P1 == g2, Class == p1) %>% nrow
#    diff = gtm - gtn
      
#    num_poss = which(int_gen$pheno == p1 &
#                      int_gen$rs_new_var1 == 5 &
#                      int_gen$rs_new_var2 == 5)
    
#    if(diff > length(num_poss)){
#      diff = length(num_poss)
#    }
#    if(length(num_poss) > 0){
#      ind = sample(num_poss, diff, replace = FALSE)
#    
#      int_gen$rs_new_var1[ind] = g1
#      int_gen$rs_new_var2[ind] = g2
#    }
    
    gtn = int_gen %>% filter(rs_new_var1 == g1 & rs_new_var2 == g2 & pheno == p1) %>% nrow
    diff = gtm - gtn

#print(table(paste(int_gen$rs_new_var1,int_gen$rs_new_var2,int_gen$pheno, sep = "_")))
#print(table(paste(mod_dat$M0P0,mod_dat$M0P1, mod_dat$Class, sep = "_")))
    print(7)
    if(diff > 0){
      if(g1 > 0){
        selec1 = c(1,2)
      } else {
        selec1 = c(0,1)
      }
      if(g2 > 0){
        selec2 = c(1,2)
      } else {
        selec2 = c(0,1)
      }
      num_poss = which(int_gen$pheno == p1 &
                                  (is.na(int_gen$rs_new_var1) | int_gen$rs_new_var1 == 5) &
                                  (is.na(int_gen$rs_new_var2) | int_gen$rs_new_var2 == 5) & 
                                  int_gen$s1 %in% selec1 &
                                  int_gen$s2 %in% selec2)

print(table(paste(int_gen$rs_new_var1,int_gen$rs_new_var2,int_gen$pheno, sep = "_")))
print(table(paste(mod_dat$M0P0,mod_dat$M0P1, mod_dat$Class, sep = "_")))
      print(8)
#      if(diff > length(num_poss)){
#        diff = length(num_poss)
#      }
#      if(length(num_poss) > 0){
#        ind = sample(num_poss, diff, replace = FALSE)
#                                  
#        int_gen$rs_new_var1[ind] = g1
#        int_gen$rs_new_var2[ind] = g2
#      }
      
      gtn = int_gen %>% filter(rs_new_var1 == g1 & rs_new_var2 == g2 & pheno == p1) %>% nrow
      diff = gtm - gtn

print(table(paste(int_gen$rs_new_var1,int_gen$rs_new_var2,int_gen$pheno, sep = "_")))
print(table(paste(mod_dat$M0P0,mod_dat$M0P1, mod_dat$Class, sep = "_")))
      print(9)
      if(diff > 0){
        num_poss = which(int_gen$pheno == p1 &
                                  (is.na(int_gen$rs_new_var1) | int_gen$rs_new_var1 == 5) &
                                  (is.na(int_gen$rs_new_var2) | int_gen$rs_new_var2 == 5))
        ind = sample(num_poss, diff, replace = FALSE)
        
        int_gen$rs_new_var1[ind] = g1
        int_gen$rs_new_var2[ind] = g2
        
      }
    }
  } else {
      print(10)
    num_poss = which(int_gen$s1 == g1 & 
                       int_gen$s2 == g2  & 
                       int_gen$pheno == p1 & 
                       (is.na(int_gen$rs_new_var1) |
                        int_gen$rs_new_var1 == 5)  & 
                       (is.na(int_gen$rs_new_var2) |
                          int_gen$rs_new_var2 == 5) )
    
    ind = sample(num_poss, gtm, replace = FALSE)
    rnge = num_poss[!num_poss %in% ind]
    int_gen$rs_new_var1[ind] = g1
    int_gen$rs_new_var2[ind] = g2
    int_gen$rs_new_var1[rnge] = 5
    int_gen$rs_new_var2[rnge] = 5
  }
print(table(paste(int_gen$rs_new_var1,int_gen$rs_new_var2,int_gen$pheno, sep = "_")))
print(table(paste(mod_dat$M0P0,mod_dat$M0P1, mod_dat$Class, sep = "_")))
 
}
print(11)
id_ind = which(bim$id %in% c(s1,s2)) %>% sort
if (id_ind[1] == which(bim$id == s1)){
    id_order = c("FID","IID","PID","MID","Sex","PHENOTYPE",
                 bim$id[1:id_ind[1]], 
                 paste0("rs_int1_",s1),
                 bim$id[(id_ind[1]+2):id_ind[2]],
                 paste0("rs_int2_",s2),
                 bim$id[(id_ind[2]+2):length(bim$id)])
} else {
    id_order = c("FID","IID","PID","MID","Sex","PHENOTYPE",
                 bim$id[1:id_ind[1]], 
                 paste0("rs_int2_",s2),
                 bim$id[(id_ind[1]+2):id_ind[2]],
                 paste0("rs_int1_",s1),
                 bim$id[(id_ind[2]+2):length(bim$id)])
}
print(12)
int_gen = int_gen %>%
    select(-pheno, -s1, -s2)
#print(head(int_gen))
print(12.1)
gen = gen %>%
    select(-rs_new_var1, -rs_new_var2) %>%
    left_join(int_gen, by = "id") %>%
    mutate(FID = id,
           IID = id,
           PID = 0,
           MID = 0,
           Sex = 1,
           PHENOTYPE = pheno)
print(12.2)
colnames(gen)[colnames(gen) == "s1"] = s1
colnames(gen)[colnames(gen) == "s2"] = s2
colnames(gen)[colnames(gen) == "rs_new_var1"] = paste0("rs_int1_",s1)
colnames(gen)[colnames(gen) == "rs_new_var2"] = paste0("rs_int2_",s2)
print(12.3)
#print(colnames(gen))
#print(id_order)

gen = gen %>%
    select(all_of(id_order)) %>%
    as.data.frame
print(13)

#for (snp in colnames(gen)[grepl("^rs",colnames(gen))]){
#    print(snp)
#    if (length(bim$id[bim$id == snp]) == 1){
#        a1 = bim$a1[bim$id == snp]
#        a2 = bim$a2[bim$id == snp]
#    } else {
#        a1 = "A"
#        a2 = "T"
#    }
#    gen_rep = as.character(unlist(gen[snp]))
#    gen_rep = gsub("2",paste0(a1,"_",a1),gen_rep)
#    gen_rep = gsub("1",paste0(a1,"_",a2),gen_rep)
#    gen_rep = gsub("0",paste0(a2,"_",a2),gen_rep)
#    gen[snp] = gen_rep

#    gen = gen %>%
#        separate(!!snp, into = paste0(snp, c(1,2)), sep = "_")
#}
print("bim")

bim = data.frame(id = colnames(gen)[grepl("^rs",colnames(gen))]) %>%
    left_join(bim, by = "id") %>%
    select(chr, id, cm, bp, a1, a2)

rownames(bim) = bim$id

print("repl")
repl = which(is.na(bim$chr))
for (r in repl){
    print(r)
    bim$chr[r] = bim$chr[r-1]
    bim$cm[r] = 0
    bim$bp[r] = bim$bp[r-1] + 1
    bim$a1[r] = bim$a1[r-1]
    bim$a2[r] = bim$a2[r-1]
}

print("snpMatrix")

snp_transform = function(snp){
  snp = as.numeric(snp)
  snp[snp == 2] = 3
  snp[snp == 1] = 2
  snp[snp == 0] = 1
  snp[is.na(snp)] = 0

  return(snp)
}

snpmat = gen %>%
    select(starts_with("rs")) 

snpmat[] = lapply(snpmat, snp_transform)

rownames(snpmat) = gen$IID

snpmat = as.matrix(snpmat) 
snpmat = new("SnpMatrix", snpmat)

fam = gen %>%
    select(1:6) %>%
    mutate(PHENOTYPE = ifelse(PHENOTYPE == 1, 2, 1))
rownames(fam) = gen$IID

output_file = gsub("_models","",args)
output_file = gsub(".tsv","",output_file)
output_file = gsub(".txt","",output_file)
output_file = paste0("chr22_",output_file)

#write_tsv(gen, "test.gen")
#write_tsv(bim, "test.bim")
#write_tsv(fam, "test.fam")

#colnames(fam) = c("pedigree", "member", "father", "mother", "sex",
#                      "affected")
#colnames(bim) = c("chromosome", "snp.name", "cM", "position",
#                       "allele.1", "allele.2")

print("write plink")
#output_file = "test"
write.plink(output_file, snps = snpmat, subject.data = fam, snp.data = bim)
write_tsv(fam, paste0(output_file,".fam"), col_names = FALSE)
write_tsv(bim, paste0(output_file,".bim"), col_names = FALSE)

print(14)
#write_tsv(gen, paste0(output_file,".ped"), col_names = FALSE)
#write_tsv(bim, paste0(output_file,".map"), col_names = FALSE)


