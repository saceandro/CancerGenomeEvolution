#!/usr/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N alpha_beta_same_cn_grad_10000
#$ -o grad_result/4_2_3_2_10000_10000.log
#$ -e grad_result/4_2_3_2_10000_10000.err

./grad 4 2 10000 4_2_3_2_10000_100000.read grad_result/4_2_3_2_10000_10000.u_n grad_result/4_2_3_2_10000_10000.llik grad_result/4_2_3_2_10000_10000.llik_final
