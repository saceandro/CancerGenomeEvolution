#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -e ../log_qsub/log
#$ -o ../log_qsub/log

source ~/.zshrc

num_of_split=$1

pa_old=../params_qsub/old
pa_test=../params_qsub/new
pa_best=../params_qsub/best
vf_dvf_test=../vf_qsub/new
llik=../llik_qsub/llik
# dx_dy=../dx_dy_qsub/$iteration

print "./mstep_qsub 2 0 $num_of_split $pa_old $pa_test $vf_dvf_test $llik $pa_best"
./mstep_qsub 2 0 $num_of_split $pa_old $pa_test $vf_dvf_test $llik $pa_best
