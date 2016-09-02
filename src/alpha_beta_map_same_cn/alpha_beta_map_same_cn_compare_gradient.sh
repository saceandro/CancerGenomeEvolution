#!/usr/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N alpha_beta_map_same_cn_compare_gradient

for ((j = 0; j < 5; j++)); do
    for ((i = 1; i <= 4; i++)); do
        for ((l = 0; l < 4; l++)); do
            ./alpha_beta_map_same_cn_compare_gradient 4 10 2 alpha_beta_same_cn_generate.random.reads10 compare_gradient_result2/u/4_10_2_random_reads.${j}_${i} compare_gradient_result2/beta/4_10_2_random_reads.${j}_${l} $j $i $l;
        done
    done
done
