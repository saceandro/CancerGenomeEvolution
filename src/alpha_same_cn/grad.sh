#!/usr/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N alpha_same_cn_grad

./grad 2 0 100000 grad.2_2_1_0.beta generated/2_2_1_0_10000_100000.read grad_result/2_2_1_0_10000_100000.u_beta grad_result/2_2_1_0_10000_100000.llik grad_result/2_2_1_0_10000_100000.llik_final
