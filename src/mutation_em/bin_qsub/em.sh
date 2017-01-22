#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N em_qsub
#$ -e ../log_qsub/em.err
#$ -o ../log_qsub/em.log

source ~/.zshrc
rm -f ../log_qsub/*
rm -f ../llik_qsub/llik
rm -f ../vf_qsub/*
rm -f ../du_dn_llik_qsub/*
rm -f ../params_qsub/old
rm -f ../params_qsub/new
rm -f ../params_qsub/best
rm -f ../rmsd_qsub/rmsd
cp ../../read_generation/generated/subtype2_topology0/params/0.4_0.8_4.x_y_u_n ../params_qsub/old

pa_true=../../read_generation/generated/subtype2_topology0/params/0.1_0.3_3.x_y_u_n
fileNum=`ls -l ../data/ | wc -l`
fileNum=$(($fileNum - 1))

qsub -N em_sync0 dummy.sh
qsub -N em_sync -hold_jid em_sync0 mstep.sh $fileNum $pa_true
