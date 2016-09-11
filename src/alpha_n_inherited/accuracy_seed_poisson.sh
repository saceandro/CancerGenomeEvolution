#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N accuracy_seed_poisson
# $1: #bp, $2: seed


../read_generation/generate_reads_poisson 4 1000000 1000 2 $1 $2 ../read_generation/generated/4_2_original.u_n generated_forpoisson/4_2_original_total1000000_td1000_cell1000000_${1}_${2}.reads 1> 4_2_original_total1000000_td1000_cell1000000_${1}_${2}.log_gen 2> 4_2_original_total1000000_td1000_cell1000000_${1}_${2}.err_gen

n=`wc -l generated_forpoisson/4_2_original_total1000000_td1000_cell1000000_${1}_${2}.reads`

./grad_foraccuracy_n 4 2 $n ../read_generation/generated/4_2_original.purity generated_forpoisson/4_2_original_total1000000_td1000_cell1000000_${1}_${2}.reads grad_result_forpoisson/4_2_original_total1000000_td1000_cell1000000_${1}_${2}.u_n grad_result_forpoisson/4_2_original_total1000000_td1000_cell1000000_${1}_${2}.llik grad_result_forpoisson/4_2_original_total1000000_td1000_cell1000000_${1}_${2}.llik_final > grad_result_forpoisson/4_2_original_total1000000_td1000_cell1000000_${1}_${2}.log

./calc_params 4 2 grad_result_forpoisson/4_2_original_total1000000_td1000_cell1000000_${1}_${2}.u_n grad_result_forpoisson/4_2_original_total1000000_td1000_cell1000000_${1}_${2}.t_n

./calc_accuracy 4 ../read_generation/generated/4_2_original.u_n ../read_generation/generated/4_2_original.t_n grad_result_forpoisson/4_2_original_total1000000_td1000_cell1000000_${1}_${2}.u_n grad_result_forpoisson/4_2_original_total1000000_td1000_cell1000000_${1}_${2}.t_n accuracy_result_seed_poisson/4_2_original_total1000000_td1000_cell1000000_${1}_${2}.rmsd.u accuracy_result_seed_poisson/4_2_original_total1000000_td1000_cell1000000_${1}_${2}.rmsd.t accuracy_result_seed_poisson/4_2_original_total1000000_td1000_cell1000000_${1}_${2}.rmsd.n
