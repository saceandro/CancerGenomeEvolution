#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -e ../log_qsub_coverage/log
#$ -o ../log_qsub_coverage/log

source ~/.zshrc

coverage=$1
num_of_split=$2
pa_true=$3

pa_old=../params_qsub_coverage/old${coverage}
pa_test=../params_qsub_coverage/new${coverage}
pa_best=../params_qsub_coverage/best${coverage}
pa_log=../params_qsub_coverage/log${coverage}
vf_dvf_test=../vf_qsub_coverage/new${coverage}
vf_dvf_old=../vf_qsub_coverage/old${coverage}
llik=../llik_qsub_coverage/llik${coverage}
rmsd=../rmsd_qsub_coverage/rmsd${coverage}
params_diff=../padiff_qsub_coverage/padiff${coverage}
# dx_dy=../dx_dy_qsub_coverage/$iteration

./mstep_qsub_coverage 2 0 $num_of_split $pa_old $pa_test $vf_dvf_test $llik $pa_best $vf_dvf_old $pa_log $pa_true $rmsd $params_diff $coverage
