#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -e ../log_qsub_coverage_seed_rand_param/esteperr
#$ -o ../log_qsub_coverage_seed_rand_param/esteplog
#$ -l s_vmem=8G,mem_req=8G
# $1: snvs
source ~/.zshrc

u1=$1
u2=$2
n=$3

snvs=$4
coverage=$5
seed=$6
iter=$7

datafile_id=$8

filename=../data_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${datafile_id}
pa_old=../params_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}/old
pa_test=../params_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}/new
vf_dvf_old=../vf_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}/old
vf_dvf_test=../vf_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}/new
du_dn_llik=../du_dn_llik_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}/$datafile_id
# errf=../log_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}/$datafile_id

filelen=`wc -l $filename`
../bin/estep_mem_multinomial 2 0 $filelen $pa_old $pa_test $vf_dvf_old $vf_dvf_test $filename $du_dn_llik # 1> ../log/e${iteration}_${SGE_TASK_ID}.log 2> ../log/e${iteration}_${SGE_TASK_ID}.err
