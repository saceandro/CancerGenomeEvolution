#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N em_qsub_snv
#$ -e ../log_qsub_snv/em.err
#$ -o ../log_qsub_snv/em.log

snvs=$1

source ~/.zshrc
cp ../../read_generation/generated/subtype2_topology0/params/0.4_0.8_4.x_y_u_n ../params_qsub_snv/old${snvs}

pa_true=../../read_generation/generated/subtype2_topology0/params/0.1_0.3_3.x_y_u_n


fileNum=`ls -l ../data/ | wc -l`
fileNum=$(($fileNum - 1))

qsub -N em_snv${snvs}_init dummy_snv.sh $snvs
qsub -N em_snv${snvs} -hold_jid em_snv${snvs}_init mstep_snv.sh $snvs $fileNum $pa_true
