#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -e ../log_qsub/log
#$ -o ../log_qsub/log

source ~/.zshrc

num_of_split=$1
pa_true=$2

pa_old=../params_qsub/old
pa_test=../params_qsub/new
pa_best=../params_qsub/best
pa_log=../params_qsub/log
vf_dvf_test=../vf_qsub/new
vf_dvf_old=../vf_qsub/old
llik=../llik_qsub/llik
rmsd=../rmsd_qsub/rmsd
params_diff=../padiff_qsub/padiff
# dx_dy=../dx_dy_qsub/$iteration

./mstep_qsub 2 0 $num_of_split $pa_old $pa_test $vf_dvf_test $llik $pa_best $vf_dvf_old $pa_log $pa_true $rmsd $params_diff
