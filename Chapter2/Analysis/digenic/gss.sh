file_in=$1

file_in=$( echo "'"${file_in}"'" )

matlab -r "gwis-stats/runGssAnalysis(${file_in})"
