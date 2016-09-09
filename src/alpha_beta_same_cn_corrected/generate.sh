#!/usr/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N generate

./generate_params 4 2 10000 3 generated_foraccuracy/4_2_3.purity generated_foraccuracy/4_2_3.u_n_xi generated_foraccuracy/4_2_3.t_n 2
# for ((n = 1; n <= 100000; n*=2)); do
#     ./generate_reads 4 2 10000 $n $n generated_foraccuracy/4_2_3.u_n_xi generated_foraccuracy/4_2_3_10000.read$n 2;
# done
