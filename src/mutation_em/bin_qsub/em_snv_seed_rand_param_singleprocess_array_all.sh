#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N em_qsub_snv_seed_rand_param_all
#$ -e ../log_qsub_snv_seed_rand_param/em.err
#$ -o ../log_qsub_snv_seed_rand_param/em.log
#$ -l s_vmem=1G,mem_req=1G

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
u2=0.4
n=1

snvs=3000
./em_snv_seed_rand_param_singleprocess_array.sh $u_lower $em_max_iter $grad_desc_max_iter $snvs $u1 $u2 $n
echo "========================================================================="

# for ((snvs=10; snvs<=100; snvs+=10)); do
#         qsub -N em_snv${snvs}_seed_rand_all -l ljob em_snv_seed_rand_init.sh $u_lower $em_max_iter $grad_desc_max_iter $snvs $iter
# done
