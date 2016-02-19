#!/bin/zsh
#$ -S /bin/zsh
#$ -cwd
#$ -N mixture_cn_accuracy
#$ -o ./mixture_cn_io
#$ -e ./mixture_cn_io

source ~/.zshrc
echo SGE_TASK_ID:$SGE_TASK_ID
echo SGE_TASK_FIRST:$SGE_TASK_FIRST
echo SGE_TASK_LAST:$SGE_TASK_LAST
echo SGE_TASK_STEPSIZE:$SGE_TASK_STEPSIZE

n=$(2 ** $SGE_TASK_ID)
./mixture_cn_generate `printf "./mixture_cn_io/mixture_cn_generate.out%d" $n` 100000 $n 4 kappa.txt
./mixture_cn_accuracy 4 4 $n `printf "./mixture_cn_io/mixture_cn_generate.out%d" $n` `printf "./mixture_cn_io/mixture_cn_accuracy.out%d" $n` kappa.txt `printf "./mixture_cn_io/mixture_cn_accuracy.accuracy%d" $n`
