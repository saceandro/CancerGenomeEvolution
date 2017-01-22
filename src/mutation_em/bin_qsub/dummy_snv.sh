#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -e ../log_qsub_snv/em0.err
#$ -o ../log_qsub_snv/em0.log
# $1: snvs

snvs=$1
param0=../params_qsub_snv/old${snvs}
vf0=../vf_qsub_snv/old${snvs}

../bin/calc_vf_dvf_multinomial 2 0 $param0 $vf0
