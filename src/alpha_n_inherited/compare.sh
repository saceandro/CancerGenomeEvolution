#!/usr/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N compare

for ((uindex = 1; uindex<=4; uindex++)); do
    for ((nindex = 1; nindex<=4; nindex++)); do
        ./compare 4 1000 2 0.97547454735 generated/4_2_3_10000_100000.reads compare_result4/u.${uindex} compare_result4/n.${nindex} 2 $uindex $nindex > compare_result4/log.${uindex}_${nindex}
    done
done
