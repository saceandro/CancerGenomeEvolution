#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N em_qsub_snv_seed_rand_param_all
#$ -e ../log_qsub_snv_seed_rand_param/em.err
#$ -o ../log_qsub_snv_seed_rand_param/em.log

source ~/.zshrc

u_lower=$1
em_max_iter=$2
grad_desc_max_iter=$3

mkdir -p ../log_qsub_snv_seed_rand_param
rm -f ../log_qsub_snv_seed_rand_param/em.log
rm -f ../log_qsub_snv_seed_rand_param/em.err
rm -f ../log_qsub_snv_seed_rand_param/esteplog
rm -f ../log_qsub_snv_seed_rand_param/esteperr

u1=0.5
u2=0.9
n=3

snvs=3000
for ((seed=1; seed<=10; ++seed)); do
        qsub -N em_snv${snvs}_seed_rand_param${u1}_${u2}_${n} -l ljob em_snv_seed_rand_param.sh $u_lower $em_max_iter $grad_desc_max_iter $snvs $seed $u1 $u2 $n
done
echo "========================================================================="

# for ((snvs=10; snvs<=100; snvs+=10)); do
#         qsub -N em_snv${snvs}_seed_rand_all -l ljob em_snv_seed_rand_init.sh $u_lower $em_max_iter $grad_desc_max_iter $snvs $iter
# done
