#!/usr/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N alpha_same_cn_compare

for ((j = 0; j < 5; j++)); do
    for ((i = 1; i <= 4; i++)); do
        ./compare 4 10 2 $j $i compare.beta alpha_beta_same_cn_generate.random.reads10 compare_result/4_10_2_random_reads.${j}_${i}
    done
done
