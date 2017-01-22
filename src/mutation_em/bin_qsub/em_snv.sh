#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N em_qsub_snv
#$ -e ../log_qsub_snv/em.err
#$ -o ../log_qsub_snv/em.log

source ~/.zshrc
rm -f ../log_qsub_snv/*
rm -f ../llik_qsub_snv/llik*
rm -f ../vf_qsub_snv/*
rm -f ../du_dn_llik_qsub_snv/*
rm -f ../params_qsub_snv/old*
rm -f ../params_qsub_snv/new*
rm -f ../params_qsub_snv/best*
rm -f ../rmsd_qsub_snv/rmsd*
rm -f ../padiff_qsub_snv/padiff*
rm -f ../data/* # remove data files

mkdir -p ../log_qsub_snv
mkdir -p ../llik_qsub_snv
mkdir -p ../vf_qsub_snv
mkdir -p ../du_dn_llik_qsub_snv
mkdir -p ../params_qsub_snv
mkdir -p ../rmsd_qsub_snv
mkdir -p ../padiff_qsub_snv

for ((snvs=10; snvs<=10000; snvs*=10)); do
    param0=../params_qsub_snv/old${snvs}
    cp ../../read_generation/generated/subtype2_topology0_snv/params/0.4_0.8_4.x_y_u_n $param0
    
    pa_true=../../read_generation/generated/subtype2_topology0_snv/params/0.1_0.3_3.x_y_u_n

    ./split.py ../../read_generation/generated/subtype2_topology0_snv/reads/0.1/0.3/3/coverage100_snv${snvs}_seed1.reads $snvs 100
    
    fileNum=`ls -l ../data/ | wc -l`
    fileNum=$(($fileNum - 1))

    vf0=../vf_qsub_snv/old${snvs}
    ../bin/calc_vf_dvf_multinomial 2 0 $param0 $vf0

    qsub -N em_snv${snvs} -sync y mstep_snv.sh $snvs $fileNum $pa_true
done
