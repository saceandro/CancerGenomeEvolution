#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -e ../log_qsub_snv_seed/esteperr
#$ -o ../log_qsub_snv_seed/esteplog
#$ -l s_vmem=8G,mem_req=8G
# $1: snvs
source ~/.zshrc

snvs=$1
seed=$2

datafile_id=$(($SGE_TASK_ID - 1))
filename=../data_snv_seed/${snvs}/${seed}/${datafile_id}
pa_old=../params_qsub_snv_seed/${snvs}/${seed}/old
pa_test=../params_qsub_snv_seed/${snvs}/${seed}/new
vf_dvf_old=../vf_qsub_snv_seed/${snvs}/${seed}/old
vf_dvf_test=../vf_qsub_snv_seed/${snvs}/${seed}/new
du_dn_llik=../du_dn_llik_qsub_snv_seed/${snvs}/${seed}/$datafile_id

n=`wc -l $filename`
../bin/estep_mem_multinomial 2 0 $n $pa_old $pa_test $vf_dvf_old $vf_dvf_test $filename $du_dn_llik # 1> ../log/e${iteration}_${SGE_TASK_ID}.log 2> ../log/e${iteration}_${SGE_TASK_ID}.err
