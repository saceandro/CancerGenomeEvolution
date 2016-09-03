#!/usr/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N accuracy

./grad 4 2 100000 grad.4_2_3_2.beta generated/4_2_3_2_10000_100000.read accuracy_result/4_2_3_2_10000_100000.u_beta accuracy_result/4_2_3_2_10000_100000.llik accuracy_result/4_2_3_2_10000_100000.llik_final accuracy.4_2_3_2.u accuracy_result/4_2_3_2_10000_100000.accuracy
