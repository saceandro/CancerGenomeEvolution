#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N accuracy_seed_bfgs2
# $1: #locus, $2: seed

./accuracy_bfgs2 4 2 $1 ../alpha_same_cn/generated_foraccuracy/4_2_3_10000.read${1}_${2} grad_result_foraccuracy_bfgs2/4_2_3_2_10000_${1}_${2}.u_n grad_result_foraccuracy_bfgs2/4_2_3_2_10000_${1}_${2}.llik grad_result_foraccuracy_bfgs2/4_2_3_2_10000_${1}_${2}.llik_final accuracy.4_2_3_2.u_n accuracy_result_seed_bfgs2/4_2_3_2_10000_${1}_${2}
