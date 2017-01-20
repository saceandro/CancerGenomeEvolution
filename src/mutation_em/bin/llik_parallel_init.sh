#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -e ../log/llik_parallel_init.err
#$ -o ../log/llik_parallel_init.log
# $1: params

params=$1
param0=../params/0
param1=../params/1
vf0=../vf/0
vf1=../vf/1

cp $params $param0
cp $params $param1
./calc_vf_dvf 2 0 $param0 $vf0
cp $vf0 $vf1
