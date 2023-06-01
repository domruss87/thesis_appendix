setwd("/Users/Dom/Documents/Bioinformatics/Epistasis/paper2/impure_models/")
set.seed(1234)

library(tidyverse)
#returns numbers of each genotype in perfect hwe distribution
hwe = function(maf, inds){
  maj_hom = (1-maf)^2
  het = 2 * maf * (1-maf)
  min_hom = maf^2
  allele_list = as.integer(as.character(c(maj_hom,het,min_hom) * inds))
  return(allele_list)
}

gt_assigner = function(s1,s2){
  jnt = paste0(s1,s2)
  tab_out = data.frame(s1AA = c(0,0,0),s1Aa = c(0,0,0),s1aa = c(0,0,0))
  tab_out[1,1] = length(jnt[jnt=="00"])
  tab_out[2,1] = length(jnt[jnt=="01"])
  tab_out[3,1] = length(jnt[jnt=="02"])
  tab_out[1,2] = length(jnt[jnt=="10"])
  tab_out[2,2] = length(jnt[jnt=="11"])
  tab_out[3,2] = length(jnt[jnt=="12"])
  tab_out[1,3] = length(jnt[jnt=="20"])
  tab_out[2,3] = length(jnt[jnt=="21"])
  tab_out[3,3] = length(jnt[jnt=="22"])
  return(tab_out)
}

#gt_tab indicates genotypes that have a causative effect, starting off with none and can be edited after
gt_tab = data.frame(s1AA = c(0,0,0),s1Aa = c(0,0,0),s1aa = c(0,0,0))
row.names(gt_tab) = c("s2BB","s2Bb","s2bb")
#sampling function
cc_sample = function(hwe_tab, gt_tab, case_with_int, inds_with_int){
  cont_with_int = inds_with_int - case_with_int
  case_wo_int = 50000 - case_with_int
  cont_wo_int = 50000 - cont_with_int
  prob_per_gt = 2 * unlist(hwe_tab) / sum(hwe_tab)
  probs_int = prob_per_gt * unlist(gt_tab)
  probs_wo_int = prob_per_gt * (1-unlist(gt_tab))
  samp_list = 1:9
  gt_case_int = sample(x = samp_list, size = case_with_int, prob = probs_int, replace = TRUE)
  gt_case_wo_int = sample(x = samp_list, size = case_wo_int, prob = probs_wo_int, replace = TRUE)
  gt_cases = c(gt_case_int, gt_case_wo_int)
  gt_cont_int = sample(x = samp_list, size = cont_with_int, prob = probs_int, replace = TRUE)
  gt_cont_wo_int = sample(x = samp_list, size = cont_wo_int, prob = probs_wo_int, replace = TRUE)
  gt_cont = c(gt_cont_int, gt_cont_wo_int)
  case_tab = data.frame(s1AA = c(length(gt_cases[gt_cases == 1]),length(gt_cases[gt_cases == 2]),length(gt_cases[gt_cases == 3])),
                        s1Aa = c(length(gt_cases[gt_cases == 4]),length(gt_cases[gt_cases == 5]),length(gt_cases[gt_cases == 6])),
                        s1aa = c(length(gt_cases[gt_cases == 7]),length(gt_cases[gt_cases == 8]),length(gt_cases[gt_cases == 9])))
  cont_tab = data.frame(s1AA = c(length(gt_cont[gt_cont == 1]),length(gt_cont[gt_cont == 2]),length(gt_cont[gt_cont == 3])),
                        s1Aa = c(length(gt_cont[gt_cont == 4]),length(gt_cont[gt_cont == 5]),length(gt_cont[gt_cont == 6])),
                        s1aa = c(length(gt_cont[gt_cont == 7]),length(gt_cont[gt_cont == 8]),length(gt_cont[gt_cont == 9])))
  generated_data = data.frame(phenotype = c(rep(1,50000), rep(0,50000)),
                              M0P0 = c(rep(0,case_tab[1,1]),rep(0,case_tab[2,1]),rep(0,case_tab[3,1]),
                                       rep(1,case_tab[1,2]),rep(1,case_tab[2,2]),rep(1,case_tab[3,2]),
                                       rep(2,case_tab[1,3]),rep(2,case_tab[2,3]),rep(2,case_tab[3,3]),
                                       rep(0,cont_tab[1,1]),rep(0,cont_tab[2,1]),rep(0,cont_tab[3,1]),
                                       rep(1,cont_tab[1,2]),rep(1,cont_tab[2,2]),rep(1,cont_tab[3,2]),
                                       rep(2,cont_tab[1,3]),rep(2,cont_tab[2,3]),rep(2,cont_tab[3,3])),
                              M0P1 = c(rep(0,case_tab[1,1]),rep(1,case_tab[2,1]),rep(2,case_tab[3,1]),
                                       rep(0,case_tab[1,2]),rep(1,case_tab[2,2]),rep(2,case_tab[3,2]),
                                       rep(0,case_tab[1,3]),rep(1,case_tab[2,3]),rep(2,case_tab[3,3]),
                                       rep(0,cont_tab[1,1]),rep(1,cont_tab[2,1]),rep(2,cont_tab[3,1]),
                                       rep(0,cont_tab[1,2]),rep(1,cont_tab[2,2]),rep(2,cont_tab[3,2]),
                                       rep(0,cont_tab[1,3]),rep(1,cont_tab[2,3]),rep(2,cont_tab[3,3]))
                              )
  return(list(generated_data,case_tab,cont_tab))
  
}

gt_tab = data.frame(s1AA = c(0,0,0),s1Aa = c(0,0,0),s1aa = c(0,0,1))
maf1 = 0.4
maf2 = 0.4
effect_size = 1.2
case_num = 50000

snp_generator = function(gt_tab,maf1,maf2,effect_size,case_num){
  hwe_s1 = hwe(maf1,2*case_num)
  hwe_s2 = hwe(maf2,2*case_num)
  s1list = c(rep(0, hwe_s1[1]), rep(1, hwe_s1[2]), rep(2,hwe_s1[3])) 
  s1list = sample(s1list, 2*case_num, replace = FALSE)
  s2list = c(rep(0, hwe_s2[1]), rep(1, hwe_s2[2]), rep(2,hwe_s2[3]))
  s2list = sample(s2list, 2*case_num, replace = FALSE)
  hwe_tab = gt_assigner(s1list,s2list)
  inds_with_int = hwe_tab[1,1] * gt_tab[1,1] + 
                      hwe_tab[1,2] * gt_tab[1,2] + 
                      hwe_tab[1,3] * gt_tab[1,3] +
                      hwe_tab[2,1] * gt_tab[2,1] +
                      hwe_tab[2,2] * gt_tab[2,2] +
                      hwe_tab[2,3] * gt_tab[2,3] +
                      hwe_tab[3,1] * gt_tab[3,1] +
                      hwe_tab[3,2] * gt_tab[3,2] +
                      hwe_tab[3,3] * gt_tab[3,3] 
  aq = effect_size - 1
  bq = inds_with_int - effect_size*case_num - effect_size*inds_with_int - case_num 
  cq = effect_size * inds_with_int * case_num
  delta = bq^2-4*aq*cq
  case_with_int = ((-bq-sqrt(delta))/(2*aq)) %>% round
  
  generated_data = cc_sample(hwe_tab, gt_tab, case_with_int, inds_with_int)
  
  return(list(generated_data,hwe_tab))      #case_tab,cont_tab,hwe_tab))
}

gt_tab = data.frame(s1AA = c(0,0,0),s1Aa = c(0,0,0),s1aa = c(0,0,1))
maf1 = 0.4
maf2 = 0.4
effect_size = 1.5
case_num = 50000

x = snp_generator(gt_tab,maf1,maf2,effect_size,case_num)

jr_aaBB_tab = data.frame(s1AA = c(0,0,0),s1Aa = c(0,0,0),s1aa = c(1,0,0))
jr_AAbb_tab = data.frame(s1AA = c(0,0,1),s1Aa = c(0,0,0),s1aa = c(0,0,0))
jr_aabb_tab = data.frame(s1AA = c(0,0,0),s1Aa = c(0,0,0),s1aa = c(0,0,1))
jr_AABB_tab = data.frame(s1AA = c(1,0,0),s1Aa = c(0,0,0),s1aa = c(0,0,0))

jd_aaBB_tab = data.frame(s1AA = c(0,0,0),s1Aa = c(1,1,0),s1aa = c(1,1,0))
jd_AAbb_tab = data.frame(s1AA = c(0,1,1),s1Aa = c(0,1,1),s1aa = c(0,0,0))
jd_aabb_tab = data.frame(s1AA = c(0,0,0),s1Aa = c(0,1,1),s1aa = c(0,1,1))
jd_AABB_tab = data.frame(s1AA = c(1,1,0),s1Aa = c(1,1,0),s1aa = c(0,0,0))

off_jr_Aabb_tab = data.frame(s1AA = c(0,0,0),s1Aa = c(0,0,1),s1aa = c(0,0,1))
off_jr_aaBb_tab = data.frame(s1AA = c(0,0,0),s1Aa = c(0,0,0),s1aa = c(1,1,0))
off_jr_AABb_tab = data.frame(s1AA = c(1,1,0),s1Aa = c(0,0,0),s1aa = c(0,0,0))
off_jr_AaBB_tab = data.frame(s1AA = c(1,0,0),s1Aa = c(1,0,0),s1aa = c(0,0,0))

off_dom_Aabb_tab = data.frame(s1AA = c(1,1,1),s1Aa = c(1,1,0),s1aa = c(0,0,0))
off_dom_AaBB_tab = data.frame(s1AA = c(0,0,0),s1Aa = c(0,1,1),s1aa = c(1,1,1))
off_dom_aaBb_tab = data.frame(s1AA = c(1,1,0),s1Aa = c(1,1,0),s1aa = c(1,0,0))
off_dom_AABb_tab = data.frame(s1AA = c(0,0,1),s1Aa = c(0,1,1),s1aa = c(0,1,1))

off_rec_aabb_tab = data.frame(s1AA = c(0,0,1),s1Aa = c(0,0,1),s1aa = c(0,1,1))
off_rec_aaBB_tab = data.frame(s1AA = c(0,0,0),s1Aa = c(1,0,0),s1aa = c(1,1,1))
off_rec_AAbb_tab = data.frame(s1AA = c(0,1,1),s1Aa = c(0,0,1),s1aa = c(0,0,1))
off_rec_AABB_tab = data.frame(s1AA = c(1,1,0),s1Aa = c(1,0,0),s1aa = c(1,0,0))

gt_tabs_all = list(jr_aaBB_tab, jr_AAbb_tab, jr_aabb_tab, jr_AABB_tab,
                   jd_aaBB_tab,jd_AAbb_tab, jd_aabb_tab, jd_AABB_tab, 
                   off_jr_Aabb_tab, off_jr_aaBb_tab, off_jr_AABb_tab, off_jr_AaBB_tab,
                   off_dom_Aabb_tab, off_dom_AaBB_tab, off_dom_aaBb_tab, off_dom_AABb_tab,
                   off_rec_aabb_tab, off_rec_aaBB_tab, off_rec_AAbb_tab, off_rec_AABB_tab)
mod_nms = c("jr1_aaBB", "jr2_AAbb", "jr3_aabb", "jr4_AABB",
            "jd1_aaBB","jd2_AAbb", "jd3_aabb", "jd4_AABB", 
            "off_jr1_Aabb", "off_jr2_aaBb", "off_jr3_AABb", "off_jr4_AaBB",
            "off_dom1_Aabb", "off_dom2_AaBB", "off_dom3_aaBb", "off_dom4_AABb",
            "off_rec1_aabb", "off_rec2_aaBB", "off_rec3_AAbb", "off_rec4_AABB")

for (t in 1:20){
  tab = gt_tabs_all[[t]]
  name_mod = mod_nms[t]
  for (es in c(1.2,1.4,1.6)){
    for (m in c(0.1,0.2,0.3,0.4)){
      for (mod in 1:10){
       df_out =  snp_generator(gt_tab = tab, maf1 = m, maf2 = m, effect_size = es, case_num = 50000)[[1]][[1]]
       f_out = paste0(name_mod,"_or",es,"_maf",m,"_it",mod,".tsv")
       print(f_out)
       write_tsv(df_out, f_out)
      }
    }
  }
}








