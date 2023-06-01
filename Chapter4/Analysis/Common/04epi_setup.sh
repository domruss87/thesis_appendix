mkdir epi

/rds/bear-apps/2019b/EL8-has/software/PLINK/1.9b_6.17-x86_64/plink \
    --bfile plink/filtered_imputed \
    --remove mother_ht_fam.tsv \
    --geno 0.05 \
    --maf 0.05 \
    --make-founders --indep-pairwise 1000 80 0.8 \
    --out epi/ld

/rds/bear-apps/2019b/EL8-has/software/PLINK/1.9b_6.17-x86_64/plink \
    --bfile plink/filtered_imputed \
    --remove mother_ht_fam.tsv \
    --extract epi/ld.prune.in \
    --geno 0.05 \
    --maf 0.05 \
    --make-bed --out epi/epi_data


