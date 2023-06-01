fin=$1
fout=$( echo $fin | rev | cut -d'/' -f1 | rev )
prefix="${2}"

mkdir ld_test


plink --bfile $1 \
    --blocks no-small-max-span \
    --blocks-max-kb 500 \
    --out ld_test/${prefix}${fout}

