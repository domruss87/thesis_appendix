for h in 
do
    for i in {1..10}
    do
        for maf in 1 2 3 4
        do
            java -jar gametes_2.2_dev.jar \
                -M " -h ${h} -a 0.${mod} -a 0.${mod} -o new_models/ord2-h-${h}-rep-${i}-maf-${mod}" -q 3 -p 1000 -t 100000000 -r ${i}
            java -jar gametes_2.2_dev.jar \
                -i new_models/ord2-h-${h}-rep-${i}-maf-${mod}*_Models.txt \
                -D "-n 0.01 -x 0.5 -a 2 -s 50000 -w 50000 -r 1 -o ord2-h-${h}-rep-${i}-maf-${mod}"
        done
    done
done
