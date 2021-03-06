#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
# $ -e ../log_qsub_snv/msteplog
# $ -o ../log_qsub_snv/msteperr

source ~/.zshrc

snvs=$1
num_of_split=$2
pa_true=$3
u_lower=$4
em_max_iter=$5
grad_desc_max_iter=$6

pa_old=../params_qsub_snv/old${snvs}
pa_test=../params_qsub_snv/new${snvs}
pa_best=../params_qsub_snv/best${snvs}
pa_log=../params_qsub_snv/log${snvs}
vf_dvf_test=../vf_qsub_snv/new${snvs}
vf_dvf_old=../vf_qsub_snv/old${snvs}
llik=../llik_qsub_snv/llik${snvs}
rmsd=../rmsd_qsub_snv/rmsd${snvs}
params_diff=../padiff_qsub_snv/padiff${snvs}
# dx_dy=../dx_dy_qsub_snv/$iteration

./mstep_qsub_snv 2 0 $num_of_split $pa_old $pa_test $vf_dvf_test $llik $pa_best $vf_dvf_old $pa_log $pa_true $rmsd $params_diff $snvs $u_lower $em_max_iter $grad_desc_max_iter 1> ../log_qsub_snv/log${snvs} 2> ../log_qsub_snv/err${snvs}
