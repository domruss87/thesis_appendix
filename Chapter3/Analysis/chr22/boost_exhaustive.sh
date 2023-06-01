i=$1
out=$(echo $i | rev | cut -d'/' -f1 | rev)

echo $i

plink --bfile $i --fast-epistasis boost --epi1 0.001 --out $out



