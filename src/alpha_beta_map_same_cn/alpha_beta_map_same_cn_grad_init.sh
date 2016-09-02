#!/usr/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N alpha_beta_map_same_cn_grad_init

for ((i = 0; i < 8; i++)); do
    ./alpha_beta_map_same_cn_grad_init 4 2 10000 generated/4_2_3_2_10000_100000.read generated/4_2_3_2.u_beta grad_init_result/alpha_beta_map_same_cn_grad_init.4_2_3_2_$i.u_beta grad_init_result/alpha_beta_map_same_cn_grad_init.4_2_3_2_$i.llik $i
done
