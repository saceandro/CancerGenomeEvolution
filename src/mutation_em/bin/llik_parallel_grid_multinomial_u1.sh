#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N llik_parallel_grid_multinomial_u1
#$ -e ../log/llik_parallel_grid_multinomial_u1.err
#$ -o ../log/llik_parallel_grid_multinomial_u1.log
# 1: iteration

source ~/.zshrc
rm -f ../log/*
rm -f ../llik/llik
rm -f ../vf/*
rm -f ../du_dn_llik/*
rm -f ../params/*
rm -f ../llik_only/llik

params_dir=../../read_generation/generated/subtype2_topology0/params/
u2=0.3
n=3

fileNum=`ls -l ../data/ | wc -l`
fileNum=$(($fileNum - 1))

qsub -N total0 llik_parallel_dummy.sh
for ((u1int=1; u1int<10; ++u1int)); do
    u1=`printf '%.1f' $((u1int * 0.1))`
    u1u2n=${u1}_${u2}_${n}
    qsub -N `printf "init%d" $u1int` -hold_jid `printf "total%d" $(($u1int - 1))` llik_parallel_init_multinomial.sh ${params_dir}/${u1u2n}.x_y_u_n
    qsub -N `printf "u%d" $u1int` -hold_jid `printf "init%d" $u1int` -tc 200 -t 1-${fileNum}:1 estep_llik_grid_multinomial.sh $u1
    qsub -N `printf "total%d" $u1int` -hold_jid `printf "u%d" $u1int` llik_total_multinomial.sh $fileNum $u1
done
