#!/bin/bash
#SBATCH --mail-type ALL
#SBATCH --nodes 1
#SBATCH --time 20:0:0
#SBATCH --mem 20G
#SBATCH --job-name trial-run
module purge; module load bluebear


plink \
    --bfile plink/filtered_imputed \
    --keep plink/cases.tsv \
    --geno 0.05 \
    --maf 0.01 \
    --make-bed --out plink/cases_only

plink --bfile plink/cases_only --exclude high-ld.txt --make-founders --indep-pairwise 1000 80 0.1 --out plink/cases_onlyLD

plink --bfile plink/cases_only --extract plink/cases_onlyLD.prune.in --nonfounders --make-bed --out plink/for_king

king -b plink/for_king.bed --related --degree 3 --prefix plink/king

flashpca --bfile plink/for_king --suffix .pca

flashpca \
    --bfile plink/for_king --check \
    --outvec plink/eigenvec.txt \
    --outval plink/eigenval.txt

rm trio_eigenvectors.pca
for i in $(grep 2$ plink/filtered_imputed.fam | cut -f2)
do
    pid=$(grep $i plink/filtered_imputed.fam | grep 2$ | cut -f3)
    mid=$(grep $i plink/filtered_imputed.fam | grep 2$ | cut -f4)

    grep $i eigenvectors.pca >> trio_eigenvectors.pca 
    grep $i eigenvectors.pca | sed "s/\t${i}\t/\t${pid}\t/g" >> trio_eigenvectors.pca 
    grep $i eigenvectors.pca | sed "s/\t${i}\t/\t${mid}\t/g" >> trio_eigenvectors.pca

done


##standard case control
/rds/bear-apps/2019b/EL8-has/software/PLINK/1.9b_6.17-x86_64/plink \
    --bfile plink/filtered_imputed \
    --covar trio_eigenvectors.pca \
    --remove mother_ht_fam.tsv \
    --geno 0.05 \
    --hwe 5e-8 \
    --maf 0.05 \
    --logistic \
    --out case_control

awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' case_control.assoc.logistic | grep ADD | grep -v NA | sort -gk9 > case_control.sort.assoc.logistic

##tdt
/rds/bear-apps/2019b/EL8-has/software/PLINK/1.9b_6.17-x86_64/plink \
    --bfile plink/filtered_imputed \
    --covar eigenvectors.pca \
    --remove mother_ht_fam.tsv \
    --geno 0.1 \
    --hwe 5e-8 \
    --maf 0.05 \
    --tdt \
    --out tdt

awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' tdt.tdt | grep -v NA | sort -gk10 > tdt.sort.tdt
