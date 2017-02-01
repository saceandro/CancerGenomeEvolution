#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -e ../log_qsub_snv_seed_rand_param/esteperr
#$ -o ../log_qsub_snv_seed_rand_param/esteplog
#$ -l s_vmem=8G,mem_req=8G
# $1: snvs
source ~/.zshrc

u1=$1
u2=$2
n=$3

snvs=$4
seed=$5
iter=$6

datafile_id=$(($SGE_TASK_ID - 1))
filename=../data_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${datafile_id}
pa_old=../params_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}/old
pa_test=../params_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}/new
vf_dvf_old=../vf_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}/old
vf_dvf_test=../vf_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}/new
du_dn_llik=../du_dn_llik_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}/$datafile_id
# errf=../log_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}/$datafile_id

filelen=`wc -l $filename`
../bin/estep_mem_multinomial 2 0 $filelen $pa_old $pa_test $vf_dvf_old $vf_dvf_test $filename $du_dn_llik # 1> ../log/e${iteration}_${SGE_TASK_ID}.log 2> ../log/e${iteration}_${SGE_TASK_ID}.err
