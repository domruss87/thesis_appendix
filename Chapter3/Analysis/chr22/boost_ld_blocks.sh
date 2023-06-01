i=$1
out=$(echo $i | rev | cut -d'/' -f1 | rev)

echo $i

../ld_test_snps.sh $i 
Rscript ../ld_test_snps.R ${out} ${i}.bim  

plink --bfile ${i} --extract ld_test/${out}_snps.txt --fast-epistasis boost --epi1 0.001 --out $out

