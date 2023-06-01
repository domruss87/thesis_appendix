t=$1
mkdir $t

plink --bfile epi_big \
        --fam ${t}.fam \
        --fast-epistasis boost --epi1 0.0005 --allow-no-sex --out ${t}/${t} 

