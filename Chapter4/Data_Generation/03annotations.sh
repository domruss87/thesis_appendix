fin=$1
fout=$(echo $fin | sed 's/.vcf.gz/_annotated.tab.gz/g')

vep --cache --dir_cache ../joint_trio_filtered/ --cache_version 106 --offline --force_overwrite --fasta /rds/projects/w/williaja-genetic-simulation/complex/bcb/genomes/Hsapiens/hg38/seq/hg38.fa --compress_output gzip -i $fin -o $fout --tab --everything --custom ../joint_trio_filtered/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz,gnomAGg,vcf,exact,0,AF,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH
