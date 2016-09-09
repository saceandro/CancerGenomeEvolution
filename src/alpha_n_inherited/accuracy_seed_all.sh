#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N accuracy_seed_all

for ((n = 1; n <= 100000; n *= 2)); do
    for ((seed = 1; seed<=10; seed++)); do
        qsub -S /usr/local/bin/zsh -cwd -N accuracy_seed_${n}_${seed} -e accuracy_log/accuracy_seed_all_${n}_${seed}.err -o accuracy_log/accuracy_seed_all_${n}_${seed}.log ./accuracy_seed.sh $n $seed
    done
done
