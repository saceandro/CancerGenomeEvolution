#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -e ../log_qsub_snv_seed_rand_param/msteplog
#$ -o ../log_qsub_snv_seed_rand_param/msteperr

source ~/.zshrc

u1=$1
u2=$2
n=$3
snvs=$4
seed=$5
iter=$6
num_of_split=$7
pa_true=$8
u_lower=$9
em_max_iter=$10
grad_desc_max_iter=$11

pa_old=../params_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}/old
pa_test=../params_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}/new
pa_best=../params_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}/best
pa_log=../params_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}/log
vf_dvf_test=../vf_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}/new
vf_dvf_old=../vf_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}/old
llik=../llik_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}
rmsd=../rmsd_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}
params_diff=../padiff_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}

./mstep_qsub_snv_seed_rand_param 2 0 $num_of_split $pa_old $pa_test $vf_dvf_test $llik $pa_best $vf_dvf_old $pa_log $pa_true $rmsd $params_diff $u1 $u2 $n $snvs $seed $iter $u_lower $em_max_iter $grad_desc_max_iter 1> ../log_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}.log 2> ../log_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}.err
