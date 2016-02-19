#!/bin/zsh
#$ -S /bin/zsh
#$ -cwd
#$ -N mixture_cn_accuracy

source ~/.zshrc
echo SGE_TASK_ID:$SGE_TASK_ID
echo SGE_TASK_FIRST:$SGE_TASK_FIRST
echo SGE_TASK_LAST:$SGE_TASK_LAST
echo SGE_TASK_STEPSIZE:$SGE_TASK_STEPSIZE

n=`expr 10 \* $SGE_TASK_ID`
./mixture_cn_generate_randomseed `printf "mixture_cn_generate.out%d" $n` 100000 $n 4 kappa.txt
./mixture_cn_accuracy 4 4 $n `printf "mixture_cn_generate.out%d" $n` `printf "mixture_cn_accuracy.out%d" $n` kappa.txt mixture_cn_accuracy.accuracy
