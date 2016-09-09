#!/usr/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N compare

for ((uindex = 1; uindex<=4; uindex++)); do
    for ((nindex = 1; nindex<=4; nindex++)); do
        ./compare 4 203 2 0.97547454735 4_2_3_10000_100_1.reads_allinherited.true_inherited compare_result3/u.${uindex} compare_result3/n.${nindex} 2 $uindex $nindex > compare_result3/log.${uindex}_${nindex}
    done
done
