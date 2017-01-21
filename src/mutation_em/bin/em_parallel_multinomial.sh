#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N em_parallel_multinomial
#$ -e ../log/em_parallel_multinomial.err
#$ -o ../log/em_parallel_multinomial.log
# 1: iteration

source ~/.zshrc
rm -f ../log/*
rm -f ../llik/llik
rm -f ../vf/*
rm -f ../du_dn_llik/*
rm -f ../dx_dy/*
rm -f ../params/*
# cp ../2_0_3_4.x_y_u_n ../params/0
# cp ../2_0_3_4.x_y_u_n ../params/1
# cp ../init_param/2_0.x_y_u_n.test3 ../params/0
# cp ../init_param/2_0.x_y_u_n.test3 ../params/1
# cp ../../read_generation/generated/subtype2_topology0/params/0.4_0.8_4.x_y_u_n ../params/0
cp ../params.0.4_0.8_4_to_0.1_0.3_3.1e-3.1e-7.iter200/72 ../params/0

iteration=$1

fileNum=`ls -l ../data/ | wc -l`
fileNum=$(($fileNum - 1))

qsub -N m0 dummy_multinomial.sh
for i in `seq 1 1 $iteration`; do
    qsub -N `printf "e%d" $i` -hold_jid `printf "m%d" $(($i - 1))` -tc 200 -t 1-${fileNum}:1 estep_multinomial.sh $i
    qsub -N `printf "m%d" $i` -hold_jid `printf "e%d" $i` mstep_multinomial.sh $i $fileNum
done
