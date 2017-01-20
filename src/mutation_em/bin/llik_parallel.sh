#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N llik_parallel
#$ -e ../log/llik_parallel.err
#$ -o ../log/llik_parallel.log
# 1: iteration

source ~/.zshrc
rm -f ../log/*
rm -f ../llik/llik
rm -f ../vf/*
rm -f ../du_dn_llik/*
rm -f ../params/*
rm -f ../llik_only/llik

genparams_dir=../../read_generation/generated
subtypes=2
topology=0
n=4

fileNum=`ls -l ../data/ | wc -l`
fileNum=$(($fileNum - 1))

qsub -N total0 llik_parallel_dummy.sh
for ((u=1; u<10; ++u)); do
    subtopoun=${subtypes}_${topology}_${u}_${n}
    qsub -N `printf "init%d" $u` -hold_jid `printf "total%d" $(($u - 1))` llik_parallel_init.sh ${genparams_dir}/${subtopoun}.x_y_u_n
    qsub -N `printf "u%d" $u` -hold_jid `printf "init%d" $u` -tc 200 -t 1-${fileNum}:1 estep_llik.sh $u
    qsub -N `printf "total%d" $u` -hold_jid `printf "u%d" $u` llik_total.sh $fileNum $u
done
