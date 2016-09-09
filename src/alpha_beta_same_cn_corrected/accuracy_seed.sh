#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N accuracy_seed
# $1: #locus, $2: seed

./generate_reads 4 2 10000 $1 $2 generated_foraccuracy/4_2_3.u_n_xi generated_foraccuracy/4_2_3_10000.read${1}_${2} 2;

./grad_foraccuracy 4 2 $1 generated_foraccuracy/4_2_3.purity generated_foraccuracy/4_2_3_10000.read${1}_${2} grad_result_foraccuracy/4_2_3_2_10000_${1}_${2}.u_n grad_result_foraccuracy/4_2_3_2_10000_${1}_${2}.llik grad_result_foraccuracy/4_2_3_2_10000_${1}_${2}.llik_final > grad_result_foraccuracy/4_2_3_2_10000_${1}_${2}.log

./calc_params 4 2 grad_result_foraccuracy/4_2_3_2_10000_${1}_${2}.u_n grad_result_foraccuracy/4_2_3_2_10000_${1}_${2}.t_n

./calc_accuracy 4 generated_foraccuracy/4_2_3.u_n_xi generated_foraccuracy/4_2_3.t_n grad_result_foraccuracy/4_2_3_2_10000_${1}_${2}.u_n grad_result_foraccuracy/4_2_3_2_10000_${1}_${2}.t_n accuracy_result_seed/4_2_3_2_10000_${1}_${2}.rmsd.u accuracy_result_seed/4_2_3_2_10000_${1}_${2}.rmsd.t accuracy_result_seed/4_2_3_2_10000_${1}_${2}.rmsd.n

# ./accuracy 4 2 $1 ../alpha_same_cn/generated_foraccuracy/4_2_3_10000.read${1}_${2} grad_result_foraccuracy/4_2_3_2_10000_${1}_${2}.u_n grad_result_foraccuracy/4_2_3_2_10000_${1}_${2}.llik grad_result_foraccuracy/4_2_3_2_10000_${1}_${2}.llik_final accuracy.4_2_3_2.u_n accuracy_result_seed/4_2_3_2_10000_${1}_${2}
