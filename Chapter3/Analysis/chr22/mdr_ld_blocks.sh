i=$1
out=$(echo $i | rev | cut -d'/' -f1 | rev)

echo $i

../ld_test_snps.sh $i 
Rscript ../ld_test_snps.R ${out} ${i}.bim  

plink --bfile ${i} --extract ld_test/${out}_snps.txt --recode A --out ${out}

Rscript make_others.R ${out}

java -jar mdr_3.0.2.jar \
    -seed=123 -min=2 -max=2 -nolandscape -discrete_significance_metric=MAXIMUM_LIKELIHOOD -saveanalysis=${out}_mdr_moore.txt ${out}.mdr.tsv

