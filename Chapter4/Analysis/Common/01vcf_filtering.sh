id=$1

mkdir temp

bcftools norm -m- ${id}_filtered.vcf.gz |\
    bcftools view -e 'FORMAT/FT!="PASS"' - |\
    bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' - |\
    bgzip -c \
    > temp/${id}_cc.vcf.gz

tabix -p vcf temp/${id}_cc.vcf.gz

bcftools view -R ../exons_genocodev39_ucsc_tableBrowser.tsv \
    temp/${id}_cc.vcf.gz \
    | bgzip -c \
    > temp/${id}_cc2.vcf.gz

mv temp/${id}_cc2.vcf.gz \
    vcf/${id}_cc.vcf.gz

rm temp/${id}_cc.vcf.gz.tbi
tabix -p vcf vcf/${id}_cc.vcf.gz

echo two!

vep --cache --cache_version 106 --offline --force_overwrite \
#    --fasta hg38.fa --compress_output gzip \
#    -i ${id}_cc.vcf.gz\
#    -o ${id}_cc_annos.vcf.gz \
#    --tab --everything \
#    --custom gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz,gnomAGg,vcf,exact,0,AF,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH


zgrep -v '##' ../${id}/${id}_cc_annos.vcf.gz | cut -f1,7,50,71 | awk '{if($3<0.01&&$4<0.01) { print $1}}' | sort -u > ../${id}/${id}_rare4cc.txt

zgrep -v '##' ../${id}/${id}_cc_annos.vcf.gz | cut -f1,7,50,71 | awk '{if($3<0.01&&$4<0.01) { print $0}}' | sort -u > ../${id}/${id}_rare_conseq_cc.txt


bcftools view --exclude ID==@../${id}/${id}_rare4pvp.txt ${id}_cc.vcf.gz \
    > vcf/${id}.vcf

