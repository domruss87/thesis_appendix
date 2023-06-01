
mod=1
for h in 0.02 0.01 0.005
do
    for i in {1..30}
    do
        java -jar gametes_2.2_dev.jar -M " -h ${h} -a 0.2 -a 0.2 -o new_models/ord2-h-${h}-maf-0.2-rep-${i}-mod-${mod}" -q 3 -p 1000 -t 100000000 -r ${i}
        java -jar gametes_2.2_dev.jar -i new_models/ord2-h-${h}-maf-0.2-rep-${i}-mod-1_Models.txt\
        -D "-n 0.01 -x 0.5 -a 2500 -s 1000 -w 1000 -r 1 -o ord2-h-${h}-maf-0.2-rep-${i}"
    done
done

