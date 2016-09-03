#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N accuracy_array

source ~/.zshrc
echo SGE_TASK_ID:$SGE_TASK_ID
echo SGE_TASK_FIRST:$SGE_TASK_FIRST
echo SGE_TASK_LAST:$SGE_TASK_LAST
echo SGE_TASK_STEPSIZE:$SGE_TASK_STEPSIZE

n=$((2 ** $SGE_TASK_ID))
./accuracy 4 2 $n grad.4_2_3_2.beta generated/4_2_3_2_10000_100000.read grad_result/4_2_3_2_10000_$n.u_beta grad_result/4_2_3_2_10000_$n.llik grad_result/4_2_3_2_10000_$n.llik_final accuracy.4_2_3_2.u accuracy_result/4_2_3_2_10000_$n
