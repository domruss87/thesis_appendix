file_in=$1
file_out=$2

cassi -mem2 -lr -lr-th 0.05 -o ${file_out}_cassi.txt ${file_in}
