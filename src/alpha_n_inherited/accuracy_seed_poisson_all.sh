#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N accuracy_seed_poisson_all

for ((bp = 30000000000; bp >= 29296875; bp /= 2)); do
    for ((seed = 1; seed<=10; seed++)); do
        qsub -S /usr/local/bin/zsh -l s_vmem=16G,mem_req=16G -cwd -N accuracy_seed_poisson_${bp}_${seed} -e accuracy_log_poisson/accuracy_seed_poisson_all_${bp}_${seed}.err -o accuracy_log_poisson/accuracy_seed_poisson_all_${bp}_${seed}.log ./accuracy_seed_poisson.sh $bp $seed
    done
done
