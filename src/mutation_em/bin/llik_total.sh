#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -e ../log/llik_total.err
#$ -o ../log/llik_total.log
# 1: iteration

source ~/.zshrc

num_of_split=$1
u=$2
llik=../llik_only/llik

./grad_desc_llik $num_of_split $u $llik
