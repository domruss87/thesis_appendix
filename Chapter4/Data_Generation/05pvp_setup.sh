id=$1

bcftools norm -m- ${id}/${id}_filtered_mendelian.vcf.gz |\
    bcftools view -s ${id} - |\
    bcftools view -e 'GT="0|0" | ALT="*" | FORMAT/FT!="PASS"' - |\
    bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' - |\
    bgzip -c \
    > ${id}/${id}_pvp.vcf.gz

tabix -p vcf ${id}/${id}_pvp.vcf.gz

bcftools view -R exons_genocodev39_ucsc_tableBrowser.tsv \
    ${id}/${id}_pvp.vcf.gz \
    | bgzip -c \
    > ${id}/${id}_pvp2.vcf.gz

mv ${id}/${id}_pvp2.vcf.gz \
    ${id}/${id}_pvp.vcf.gz

rm ${id}/${id}_pvp.vcf.gz.tbi
tabix -p vcf ${id}/${id}_pvp.vcf.gz

vep --cache --cache_version 106 --offline --force_overwrite \
    --fasta hg38.fa --compress_output gzip \
    -i ${id}/${id}_pvp.vcf.gz\
    -o ${id}/${id}_annos.vcf.gz \
    --tab --everything \
    --custom gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz,gnomAGg,vcf,exact,0,AF,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH


zgrep -v '##' ${id}/${id}_annos.vcf.gz | cut -f1,7,51,71 | awk '{if($3<0.01&&$4<0.01) { print $1}}' | sort -u > ${id}/${id}_rare4pvp.txt

zgrep -v '##' ${id}/${id}_annos.vcf.gz | cut -f1,7,51,71 | awk '{if($3<0.01&&$4<0.01) { print $0}}' | sort -u > ${id}/${id}_rare_conseq.txt


bcftools view --include ID==@${id}/${id}_rare4pvp.txt ${id}/${id}_pvp.vcf.gz \
    > pvp/input/${id}.vcf

