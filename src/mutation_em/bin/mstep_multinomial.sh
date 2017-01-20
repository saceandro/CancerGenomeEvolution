#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -e ../log/log
#$ -o ../log/log
# 1: iteration

source ~/.zshrc

iteration=$1
num_of_split=$2

datafile_id=$(($SGE_TASK_ID - 1))
pa_old=../params/$iteration
pa_test=../params/$(($iteration + 1))
vf_dvf_test=../vf/$(($iteration + 1))
llik=../llik/llik
dx_dy=../dx_dy/$iteration

# pa_old=`printf "../params/%d" $iteration`
# pa_test=`printf "../params/%d" $(($iteration + 1))`
# vf_dvf_test=`printf "../vf/%d" $(($iteration + 1))`
# llik=`printf "../llik/llik"`

./grad_desc_multinomial 2 0 $num_of_split $pa_old $pa_test $vf_dvf_test $llik $dx_dy 1e-4 1e-7 # 1> ../log/m${iteration}_${SGE_TASK_ID}.log 2> ../log/m${iteration}_${SGE_TASK_ID}.err
