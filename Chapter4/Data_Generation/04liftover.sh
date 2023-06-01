
id=$1
fout=$(echo $id | sed 's/.vcf/_b37.vcf/g')
rej=$(echo $id | sed 's/.vcf/_rej.vcf/g')

picard LiftoverVcf \
     I=${id} \
     O=${fout} \
     CHAIN=hg38ToHg19.over.chain.gz \
     REJECT=${rej} \
     R=hg19.fa.gz
