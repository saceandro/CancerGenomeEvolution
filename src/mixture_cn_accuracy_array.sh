#!/usr/local/bin/zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N mixture_cn_accuracy
#$ -l s_vmem=8G,mem_req=8G

echo SGE_TASK_ID:$SGE_TASK_ID
echo SGE_TASK_FIRST:$SGE_TASK_FIRST
echo SGE_TASK_LAST:$SGE_TASK_LAST
echo SGE_TASK_STEPSIZE:$SGE_TASK_STEPSIZE

n=`expr 10 \* $SGE_TASK_ID`
./mixture_cn_generate `printf "mixture_cn_generate.out%d" $a` 100000 $n 4 kappa.txt
./mixture_cn_accuracy 4 4 $n `printf "mixture_cn_generate.out%d" $a` `printf "mixture_cn_accuracy.out%d" $a` kappa.txt mixture_cn_accuracy.accuracy
