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

params_dir=../../read_generation/generated/subtype2_topology0/params/
u1=0.1
n=4

fileNum=`ls -l ../data/ | wc -l`
fileNum=$(($fileNum - 1))

qsub -N total0 llik_parallel_dummy.sh
for ((u2int=1; u2int<10; ++u2int)); do
    u2=`printf '%.1f' $((u2int * 0.1))`
    u1u2n=${u1}_${u2}_${n}
    qsub -N `printf "init%d" $u` -hold_jid `printf "total%d" $(($u - 1))` llik_parallel_init.sh ${params_dir}/${u1u2n}.x_y_u_n
    qsub -N `printf "u%d" $u` -hold_jid `printf "init%d" $u` -tc 200 -t 1-${fileNum}:1 estep_llik_grid.sh $u2
    qsub -N `printf "total%d" $u` -hold_jid `printf "u%d" $u` llik_total.sh $fileNum $u2
done
