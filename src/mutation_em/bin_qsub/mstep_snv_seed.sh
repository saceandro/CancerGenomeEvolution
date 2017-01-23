#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -e ../log_qsub_snv_seed/msteplog
#$ -o ../log_qsub_snv_seed/msteperr

source ~/.zshrc

snvs=$1
seed=$2
num_of_split=$3
pa_true=$4
u_lower=$5
em_max_iter=$6
grad_desc_max_iter=$7

pa_old=../params_qsub_snv_seed/${snvs}/${seed}/old
pa_test=../params_qsub_snv_seed/${snvs}/${seed}/new
pa_best=../params_qsub_snv_seed/${snvs}/${seed}/best
pa_log=../params_qsub_snv_seed/${snvs}/${seed}/log
vf_dvf_test=../vf_qsub_snv_seed/${snvs}/${seed}/new
vf_dvf_old=../vf_qsub_snv_seed/${snvs}/${seed}/old
llik=../llik_qsub_snv_seed/${snvs}/${seed}
rmsd=../rmsd_qsub_snv_seed/${snvs}/${seed}
params_diff=../padiff_qsub_snv_seed/${snvs}/${seed}

./mstep_qsub_snv_seed 2 0 $num_of_split $pa_old $pa_test $vf_dvf_test $llik $pa_best $vf_dvf_old $pa_log $pa_true $rmsd $params_diff $snvs $seed $u_lower $em_max_iter $grad_desc_max_iter 1> ../log_qsub_snv_seed/${snvs}/${seed}.log 2> ../log_qsub_snv_seed/${snvs}/${seed}.err
