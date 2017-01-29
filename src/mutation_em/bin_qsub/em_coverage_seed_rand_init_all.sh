#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N em_qsub_coverage_seed_rand_all
#$ -e ../log_qsub_coverage_seed_rand/em.err
#$ -o ../log_qsub_coverage_seed_rand/em.log

source ~/.zshrc

u_lower=$1
em_max_iter=$2
grad_desc_max_iter=$3
iter=1

mkdir -p ../log_qsub_coverage_seed_rand
rm -f ../log_qsub_coverage_seed_rand/em.log
rm -f ../log_qsub_coverage_seed_rand/em.err
rm -f ../log_qsub_coverage_seed_rand/esteplog
rm -f ../log_qsub_coverage_seed_rand/esteperr

for ((coverages=10; coverages<=100; coverages+=10)); do
        qsub -N em_coverage${coverages}_seed_rand_all em_coverage_seed_rand_init.sh $u_lower $em_max_iter $grad_desc_max_iter $coverages $iter
done
