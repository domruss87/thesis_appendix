file_in=$1
file_out=$2

java -jar mdr_3.0.2.jar -seed=123 -min=3 -max=3 -cv=10 -table_data=true -saveanalysis=${file_out} ${file_in}
