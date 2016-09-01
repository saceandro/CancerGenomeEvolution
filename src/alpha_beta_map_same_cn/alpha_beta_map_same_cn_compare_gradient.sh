#!/usr/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N alpha_beta_map_same_cn_compare_gradient

for ((j = 0; j < 5; j++)); do
    for ((i = 1; i <= 4; i++)); do
        ./alpha_beta_map_same_cn_compare_gradient 4 10 2 alpha_beta_same_cn_generate.random.reads10 result/alpha_beta_map_same_cn_compare_gradient.4_10_2_random_reads.${j}_${i}.u $j $i;
    done
done
