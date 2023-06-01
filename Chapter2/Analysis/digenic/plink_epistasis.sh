file_in=$1
file_out=$2

plink --bfile $file_in --epistasis --epi1 0.01 --allow-no-sex --out ${file_out}
