file_in=$1

file_out=$2


mkdir -p ${file_out}/digenic

cd ${file_out}/digenic

java -jar SNPRuler/rule.jar 1500 2 0.2 ${file_in}

cp interactions.txt ${file_out}
