#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -e ../log_qsub_coverage/log
#$ -o ../log_qsub_coverage/log
#$ -l s_vmem=8G,mem_req=8G
# $1: coverage
source ~/.zshrc

coverage=$1

datafile_id=$(($SGE_TASK_ID - 1))
filename=../data_coverage/$datafile_id
pa_old=../params_qsub_coverage/old${coverage}
pa_test=../params_qsub_coverage/new${coverage}
vf_dvf_old=../vf_qsub_coverage/old${coverage}
vf_dvf_test=../vf_qsub_coverage/new${coverage}
du_dn_llik=../du_dn_llik_qsub_coverage/$datafile_id

n=`wc -l $filename`
print "../bin/estep_mem_multinomial 2 0 $n $pa_old $pa_test $vf_dvf_old $vf_dvf_test $filename $du_dn_llik"
../bin/estep_mem_multinomial 2 0 $n $pa_old $pa_test $vf_dvf_old $vf_dvf_test $filename $du_dn_llik # 1> ../log/e${iteration}_${SGE_TASK_ID}.log 2> ../log/e${iteration}_${SGE_TASK_ID}.err
