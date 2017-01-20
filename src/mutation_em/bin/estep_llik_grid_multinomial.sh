#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
# #$ -e ../log/estep_llik.err
# #$ -o ../log/estep_llik.log
#$ -l s_vmem=16G,mem_req=16G
# 1: iteration

source ~/.zshrc

u=$1
iteration=1

datafile_id=$(($SGE_TASK_ID - 1))
filename=../data/$datafile_id
pa_old=../params/$(($iteration-1))
pa_test=../params/$iteration
vf_dvf_old=../vf/$(($iteration-1))
vf_dvf_test=../vf/$iteration
du_dn_llik=../du_dn_llik/$datafile_id

mkdir -p ../resp/u${u}
resp=../resp/u${u}/$datafile_id

# filename=`printf "../data/%d" $datafile_id`
# pa_old=`printf "../params/%d" $(($iteration-1))`
# pa_test=`printf "../params/%d" $iteration`
# vf_dvf_old=`printf "../vf/%d" $(($iteration-1))`
# vf_dvf_test=`printf "../vf/%d" $iteration`
# du_dn_llik=`printf "../du_dn_llik/%d" $datafile_id`

n=`wc -l $filename`
./estep_mem_resp_grid_multinomial 2 0 $n $pa_old $pa_test $vf_dvf_old $vf_dvf_test $filename $du_dn_llik $resp
