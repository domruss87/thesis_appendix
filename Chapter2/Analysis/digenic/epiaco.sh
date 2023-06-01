file_in=$1

file=$(echo "'"${file_in}"'")

matlab -r "epiACO(${file},20000,250,0.2)"
