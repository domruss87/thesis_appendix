i=$1

out=$(echo $i | rev | cut -d'/' -f1 | rev)

echo $i

plink --bfile $i --logistic --out $out

awk '{if($9<0.05){print $2}}' ${out}.assoc.logistic > ${out}.snps

mkdir sig not_sig

plink --bfile $i --extract ${out}.snps --make-bed --out sig/${out}
plink --bfile $i --exclude ${out}.snps --make-bed --out not_sig/${out}

../ld_test_snps.sh sig/${out} sig
Rscript ../ld_test_snps.R sig${out} sig/${out}.bim  

plink --bfile sig/${out} --extract ld_test/sig${out}_snps.txt --recode A --out sig/${out}

Rscript make_others.R sig/${out}

../ld_test_snps.sh $i not_sig
Rscript ../ld_test_snps.R not_sig_${out} not_sig/${out}.bim  

plink --bfile not_sig/${out} --extract ld_test/not_sig_${out}_snps.txt --fast-epistasis boost --epi1 0.001 --out $out

java -jar mdr_3.0.2.jar \
    -seed=123 -min=2 -max=2 -nolandscape -discrete_significance_metric=MAXIMUM_LIKELIHOOD -saveanalysis=${out}_mdr_moore.txt sig/${out}.mdr.tsv
