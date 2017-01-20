#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N llik_parallel_grid_multinomial_n
#$ -e ../log/llik_parallel_grid_multinomial_n.err
#$ -o ../log/llik_parallel_grid_multinomial_n.log
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
u2=0.3

fileNum=`ls -l ../data/ | wc -l`
fileNum=$(($fileNum - 1))

qsub -N total0 llik_parallel_dummy.sh
for ((n=1; n<10; ++n)); do
    u1u2n=${u1}_${u2}_${n}
    qsub -N `printf "init%d" $n` -hold_jid `printf "total%d" $(($n - 1))` llik_parallel_init_multinomial.sh ${params_dir}/${u1u2n}.x_y_u_n
    qsub -N `printf "u%d" $n` -hold_jid `printf "init%d" $n` -tc 200 -t 1-${fileNum}:1 estep_llik_grid_multinomial.sh $n
    qsub -N `printf "total%d" $n` -hold_jid `printf "u%d" $n` llik_total_multinomial.sh $fileNum $n
done
