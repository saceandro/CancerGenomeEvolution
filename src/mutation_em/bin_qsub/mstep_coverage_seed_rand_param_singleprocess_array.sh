#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -e ../log_qsub_coverage_seed_rand_param/msteplog
#$ -o ../log_qsub_coverage_seed_rand_param/msteperr
# #$ -l s_vmem=4G,mem_req=4G

source ~/.zshrc

u1=$1
u2=$2
n=$3
snvs=$4
taskid_1=$(($SGE_TASK_ID - 1))
sp=$(($taskid_1 / 50))
if [[ $sp -eq 0 ]]; then
    coverage=10
elif [[ $sp -eq 1 ]]; then
    coverage=50
elif [[ $sp -eq 2 ]]; then
    coverage=100
else
    coverage=200
fi
res=$(($taskid_1 - 50 * $sp))
seed=$(($res / 5 + 1))
iter=$(($res % 5 + 1))
echo "coverage: " $coverage "seed: " $seed "iter:" $iter

num_of_split=$5
pa_true=$6
u_lower=$7
em_max_iter=$8
grad_desc_max_iter=$9

pa_old=../params_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}/old
pa_test=../params_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}/new
pa_best=../params_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}/best
pa_log=../params_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}/log
vf_dvf_test=../vf_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}/new
vf_dvf_old=../vf_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}/old
llik=../llik_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}
rmsd=../rmsd_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}
params_diff=../padiff_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}

./mstep_qsub_coverage_seed_rand_param_singleprocess 2 0 $num_of_split $pa_old $pa_test $vf_dvf_test $llik $pa_best $vf_dvf_old $pa_log $pa_true $rmsd $params_diff $u1 $u2 $n $snvs $coverage $seed $iter $u_lower $em_max_iter $grad_desc_max_iter 1> ../log_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}.log 2> ../log_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}.err
