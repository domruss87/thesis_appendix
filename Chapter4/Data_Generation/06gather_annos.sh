id=$1

for i in $(zgrep "0/0+0/0" ${id}/${id}_filtered_mendelian.vcf.gz | cut -f3); 
do 
    zgrep $i ${id}/${id}_filtered_mendelian_annotated.tab.gz | awk -v v1=${id} '{print v1,$1,$7,$19,$50,$71}' | sort -u >> ${id}/${id}_candidates.txt ; 
done ; 

