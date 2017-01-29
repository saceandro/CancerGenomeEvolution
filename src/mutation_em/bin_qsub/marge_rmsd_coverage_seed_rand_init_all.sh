#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N marge_rmsd_coverage

source ~/.zshrc

rmsd_dir=../rmsd_qsub_coverage_seed_rand_marged
mkdir -p $rmsd_dir
rmsd=${rmsd_dir}/rmsd

rm -f $rmsd

for ((coverage=10; coverage<=100; coverage+=10)); do
    for ((seed=1; seed<=10; ++seed)); do
        echo -n "$coverage\t" >> ${rmsd}
        tail -n 1 ../rmsd_qsub_coverage_seed_rand/${coverage}/${seed}/1 >> ${rmsd}
    done
done
