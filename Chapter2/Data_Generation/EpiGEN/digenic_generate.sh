for i in $(ls models/ | grep -v tri)
do

python3 simulate_data.py --corpus-id 1 --pop ASW --inds 2000 --snps 2500 --num-sims 30 --noise-maf-range 0.05 0.5 --disease-maf-range 0.39 0.4 --model ${i} --seed 123

done
