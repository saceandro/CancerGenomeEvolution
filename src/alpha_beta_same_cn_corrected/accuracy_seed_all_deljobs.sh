#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N accuracy_seed_all_bfgs2

for ((n = 1; n <= 100000; n *= 2)); do
    for ((seed = 1; seed<=10; seed++)); do
        qdel accuracy_seed_bfgs2_${n}_${seed};
    done
done
