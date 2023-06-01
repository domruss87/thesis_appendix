regenie \
    --step 2 \
    --bgen ${geno_dir}/${snp_file} --with-bgi \
    --sample ${geno_dir}/${samp} \
    --phenoFile all.phen \
    --covarFile confounderFile.tsv \
    --keep ${geno_dir}/call/biobank_qcd.id \
    --bt \
    --firth 0.01 --approx --firth-se \
    --pred all_step1_pred.list \
    --lowmem tmp \
    --bsize 400 \
    --split \
    --threads 5 \
    --out results/all_chr${chr}_${pt}

