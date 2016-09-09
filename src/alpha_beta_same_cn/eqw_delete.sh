#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N restart_jobs

for line in `cat eqw_jobs`
do
    qdel $line
done
