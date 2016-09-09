#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N accuracy_seed_all_bfgs2

for ((n = 1; n <= 100000; n *= 2)); do
    for ((seed = 1; seed<=10; seed++)); do
        qsub -S /usr/local/bin/zsh -cwd -N accuracy_seed_bfgs2_${n}_${seed} -e accuracy_log_bfgs2/accuracy_seed_all_bfgs2_${n}_${seed} -o accuracy_log_bfgs2/accuracy_seed_all_bfgs2_${n}_${seed} ./accuracy_seed_bfgs2.sh $n $seed
    done
done
