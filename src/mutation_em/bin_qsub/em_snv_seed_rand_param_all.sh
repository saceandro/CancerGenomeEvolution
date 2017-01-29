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

# mkdir -p ../log_qsub_snv_seed_rand
# rm -f ../log_qsub_snv_seed_rand/em.log
# rm -f ../log_qsub_snv_seed_rand/em.err
# rm -f ../log_qsub_snv_seed_rand/esteplog
# rm -f ../log_qsub_snv_seed_rand/esteperr

u1=0.5
u2=0.1
n=0.3

snvs=3000
for ((seed=1; seed<=10; ++seed)); do
        qsub -N em_snv${snvs}_seed_rand_param_all -l ljob em_snv_seed_rand_param.sh $u_lower $em_max_iter $grad_desc_max_iter $snvs $u1 $u2 $n
done
echo "========================================================================="

# for ((snvs=10; snvs<=100; snvs+=10)); do
#         qsub -N em_snv${snvs}_seed_rand_all -l ljob em_snv_seed_rand_init.sh $u_lower $em_max_iter $grad_desc_max_iter $snvs $iter
# done
