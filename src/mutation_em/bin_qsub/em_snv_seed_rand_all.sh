#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N em_qsub_snv_seed_rand_all
#$ -e ../log_qsub_snv_seed_rand/em.err
#$ -o ../log_qsub_snv_seed_rand/em.log

source ~/.zshrc

u_lower=$1
em_max_iter=$2
grad_desc_max_iter=$3

mkdir -p ../log_qsub_snv_seed_rand
rm -f ../log_qsub_snv_seed_rand/em.log
rm -f ../log_qsub_snv_seed_rand/em.err
rm -f ../log_qsub_snv_seed_rand/esteplog
rm -f ../log_qsub_snv_seed_rand/esteperr

for ((snvs=1000; snvs<=5000; snvs+=1000)); do
        qsub -N em_snv${snvs}_seed_rand_all em_snv_seed_rand.sh $u_lower $em_max_iter $grad_desc_max_iter $snvs
done
