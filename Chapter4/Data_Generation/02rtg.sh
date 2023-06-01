id=$1

bcftools view -f PASS \
     ${id}/${id}.joint.bcf |\
    bgzip -c \
    > ${id}/${id}_filtered.vcf.gz

rtg mendelian -i ${id}/${id}_filtered.vcf.gz -o ${id}/${id}_filtered_mendelian.vcf.gz --pedigree=family_relationships.ped -t GRCh38_no_alt_analysis_set.sdf | tee input_rtg_output.txt

plink --vcf ${id}/${id}.recode.vcf --genome --out ${id}/${id}
