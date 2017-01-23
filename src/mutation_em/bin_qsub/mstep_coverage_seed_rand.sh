#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -e ../log_qsub_coverage_seed_rand/msteplog
#$ -o ../log_qsub_coverage_seed_rand/msteperr

source ~/.zshrc

coverages=$1
seed=$2
iter=$3
num_of_split=$4
pa_true=$5
u_lower=$6
em_max_iter=$7
grad_desc_max_iter=$8

pa_old=../params_qsub_coverage_seed_rand/${coverages}/${seed}/${iter}/old
pa_test=../params_qsub_coverage_seed_rand/${coverages}/${seed}/${iter}/new
pa_best=../params_qsub_coverage_seed_rand/${coverages}/${seed}/${iter}/best
pa_log=../params_qsub_coverage_seed_rand/${coverages}/${seed}/${iter}/log
vf_dvf_test=../vf_qsub_coverage_seed_rand/${coverages}/${seed}/${iter}/new
vf_dvf_old=../vf_qsub_coverage_seed_rand/${coverages}/${seed}/${iter}/old
llik=../llik_qsub_coverage_seed_rand/${coverages}/${seed}/${iter}
rmsd=../rmsd_qsub_coverage_seed_rand/${coverages}/${seed}/${iter}
params_diff=../padiff_qsub_coverage_seed_rand/${coverages}/${seed}/${iter}

./mstep_qsub_coverage_seed_rand 2 0 $num_of_split $pa_old $pa_test $vf_dvf_test $llik $pa_best $vf_dvf_old $pa_log $pa_true $rmsd $params_diff $coverages $seed $iter $u_lower $em_max_iter $grad_desc_max_iter 1> ../log_qsub_coverage_seed_rand/${coverages}/${seed}/${iter}.log 2> ../log_qsub_coverage_seed_rand/${coverages}/${seed}/${iter}.err
