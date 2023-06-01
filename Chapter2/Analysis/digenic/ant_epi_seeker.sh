file_in=$1

new_folder=$2

snps=2500
small=$( echo "$snps * 0.1" | bc )
large=$( echo "$small * 0.5" | bc )

mkdir -p params/${new_folder}/digenic

cp parameters_digenic.txt params/${new_folder}/digenic/parameters.txt

echo "INPFILE ${file_in}.aes.csv     //input file name for case-control genotype data
OUTFILE ${file_in}.aes.digenic   //output file name for detected epistatic interactions" >> params/${new_folder}/digenic/parameters.txt

cd params/${new_folder}/digenic
sed -i "s/snps_in/${snps}/g" parameters.txt
sed -i "s/small_in/${small}/g" parameters.txt
sed -i "s/large_in/${large}/g" parameters.txt

bash ../../../AntEpiSeeker
