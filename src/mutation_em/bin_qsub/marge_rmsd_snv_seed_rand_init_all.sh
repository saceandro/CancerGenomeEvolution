#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N marge_rmsd

source ~/.zshrc

rmsd_dir=../rmsd_qsub_snv_seed_rand_marged
mkdir -p $rmsd_dir
rmsd=${rmsd_dir}/rmsd

rm -f $rmsd

for ((snvs=1000; snvs<=5000; snvs+=1000)); do
    for ((seed=1; seed<=10; ++seed)); do
        echo -n "$snvs\t" >> ${rmsd}
        tail -n 1 ../rmsd_qsub_snv_seed_rand/${snvs}/${seed}/1 >> ${rmsd}
    done
done
