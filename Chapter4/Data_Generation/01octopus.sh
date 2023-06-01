id=$1

ref="hg38.fa"
bam="alignment_files"

mkdir tmp
samtools sort -@ 10 -O BAM -T tmp/${id}_sort.bam -o ${bam}/${id}.sort.aligned.bam ${bam}/${id}.aligned.bam
mv ${bam}/${id}.sort.aligned.bam ${bam}/${id}.aligned.bam

samtools sort -@ 10 -O BAM -T tmp/${pid}_sort.bam -o ${bam}/${pid}.sort.aligned.bam ${bam}/${pid}.aligned.bam
mv ${bam}/${pid}.sort.aligned.bam ${bam}/${pid}.aligned.bam

samtools sort -@ 10 -O BAM -T tmp/${mid}_sort.bam -o ${bam}/${mid}.sort.aligned.bam ${bam}/${mid}.aligned.bam
mv ${bam}/${mid}.sort.aligned.bam ${bam}/${mid}.aligned.bam

samtools index ${bam}/${id}.aligned.bam
samtools index ${bam}/${pid}.aligned.bam
samtools index ${bam}/${mid}.aligned.bam

octopus -R ${ref} -I ${bam}/${id}.aligned.bam ${bam}/${pid}.aligned.bam ${bam}/${mid}.aligned.bam --threads --refcall POSITIONAL -o ${id}/${id}.pop.bcf

