#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N grad_n
# $1: #locus, $2: seed

./grad_foraccuracy_n 4 2 2470 ../read_generation/generated/4_2_original.purity ../read_generation/generated/4_2_original_poisson_1000000_1000_1_cell1000000.reads grad_result/4_2_original_poisson_1000000_1000_1_mr_0.5_cell1000000_betatild0.50.u_n grad_result/4_2_original_poisson_1000000_1000_1_mr_0.5_cell1000000_betatild0.5.llik grad_result/4_2_original_poisson_1000000_1000_1_mr_0.5_cell1000000_betatild0.5.llik_final > grad_result/4_2_original_poisson_1000000_1000_1_mr_0.5_cell1000000_betatild0.5.log

./calc_params 4 2 grad_result/4_2_original_poisson_1000000_1000_1_mr_0.5_cell1000000_betatild0.5.u_n grad_result/4_2_original_poisson_1000000_1000_1_mr_0.5_cell1000000_betatild0.5.t_n

./calc_accuracy 4 ../read_generation/generated/4_2_original.u_n ../read_generation/generated/4_2_original.t_n grad_result/4_2_original_poisson_1000000_1000_1_mr_0.5_cell1000000_betatild0.5.u_n grad_result/4_2_original_poisson_1000000_1000_1_mr_0.5_cell1000000_betatild0.5.t_n grad_result/4_2_original_poisson_1000000_1000_1_mr_0.5_cell1000000_betatild0.5.rmsd.u grad_result/4_2_original_poisson_1000000_1000_1_mr_0.5_cell1000000_betatild0.5.rmsd.t grad_result/4_2_original_poisson_1000000_1000_1_mr_0.5_cell1000000_betatild0.5.rmsd.n
