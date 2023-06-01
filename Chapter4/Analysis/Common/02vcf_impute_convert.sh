bcftools merge vcf/*vcf.gz -Oz -o vcf/all_trios.vcf.gz

bcftools view -t chr${chr} vcf/all_trios.vcf.gz -Oz -o vcf/all_trios_chr${chr}.vcf.gz

java -jar beagle.22Jul22.46e.jar \
    gt=vcf/all_trios_chr${chr}.vcf.gz \
    ref=panel/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel${x}.vcf.gz\
    out=vcf/imputed_chr${chr}.vcf.gz \
    map=plink.chr${chr}.GRCh38_mod.map

zgrep '^#' vcf/imputed_chr1.vcf.gz.vcf.gz > vcf/merged_imputed.vcf
for i in {1..22}; do zgrep -v '^#' vcf/imputed_chr${i}.vcf.gz.vcf.gz >> vcf/merged_imputed.vcf ; done
bgzip vcf/merged_imputed.vcf
tabix -p vcf vcf/merged_imputed.vcf.gz

bcftools view -e 'INFO/DR2<=0.8' --max-alleles 2 --min-alleles 2 vcf/merged_imputed.vcf.gz | bgzip -c > vcf/filtered_imputed.vcf.gz

tabix -p vcf vcf/filtered_imputed.vcf.gz

/rds/bear-apps/2019b/EL8-has/software/PLINK/1.9b_6.17-x86_64/plink --vcf vcf/filtered_imputed.vcf.gz --make-bed --out plink/filtered_imputed


