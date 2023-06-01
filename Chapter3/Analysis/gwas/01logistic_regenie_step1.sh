regenie \
	--step 1 \
	--bed biobank \
	--extract biobank_qcd.snplist \
	--keep biobank_qcd.id \
	--phenoFile all.phen \
	--covarFile confounderFile.tsv \
	--bsize 1000 \
    --bt \
	--lowmem tmp \
	--threads 20 \
	--out all_step1 
 

