#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -e ../log_qsub/em0.err
#$ -o ../log_qsub/em0.log

param0=../params_qsub/old
vf0=../vf_qsub/old

../bin/calc_vf_dvf_multinomial 2 0 $param0 $vf0
