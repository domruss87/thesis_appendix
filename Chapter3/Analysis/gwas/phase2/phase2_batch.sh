t=$1
mkdir -p ${t}/plink_exhaustive


plink --bfile ${t}/${t} \
        --fam ${t}.fam \
        --fast-epistasis boost --epi1 0.0005 --allow-no-sex --out ${t}/plink_exhaustive/${t}

