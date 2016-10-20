#!/usr/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N compare

for ((uindex = 1; uindex<=4; uindex++)); do
    for ((nindex = 1; nindex<=4; nindex++)); do
        ./compare 4 3427 2 0.97547454735 4_2_3_manip_10000_1000_1_cell10000.reads compare_result/u.${uindex} compare_result/n.${nindex} 2 $uindex $nindex > compare_result/log.${uindex}_${nindex}
    done
done
