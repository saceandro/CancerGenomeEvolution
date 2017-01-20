#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -e ../log/m0.err
#$ -o ../log/m0.log

param0=../params/0
param1=../params/1
vf0=../vf/0
vf1=../vf/1

cp $param0 $param1
./calc_vf_dvf_multinomial 2 0 $param0 $vf0
cp $vf0 $vf1
