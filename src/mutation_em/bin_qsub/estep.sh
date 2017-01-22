#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -e ../log_qsub/log
#$ -o ../log_qsub/log
#$ -l s_vmem=8G,mem_req=8G

source ~/.zshrc

datafile_id=$(($SGE_TASK_ID - 1))
filename=../data/$datafile_id
pa_old=../params_qsub/old
pa_test=../params_qsub/new
vf_dvf_old=../vf_qsub/old
vf_dvf_test=../vf_qsub/new
du_dn_llik=../du_dn_llik_qsub/$datafile_id

n=`wc -l $filename`
print "../bin/estep_mem_multinomial 2 0 $n $pa_old $pa_test $vf_dvf_old $vf_dvf_test $filename $du_dn_llik"
../bin/estep_mem_multinomial 2 0 $n $pa_old $pa_test $vf_dvf_old $vf_dvf_test $filename $du_dn_llik # 1> ../log/e${iteration}_${SGE_TASK_ID}.log 2> ../log/e${iteration}_${SGE_TASK_ID}.err
