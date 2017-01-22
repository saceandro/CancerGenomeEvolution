#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N em_qsub_coverage
#$ -e ../log_qsub_coverage/em.err
#$ -o ../log_qsub_coverage/em.log

source ~/.zshrc
rm -f ../log_qsub_coverage/*

mkdir -p ../log_qsub_coverage
mkdir -p ../llik_qsub_coverage
mkdir -p ../vf_qsub_coverage
mkdir -p ../du_dn_llik_qsub_coverage
mkdir -p ../params_qsub_coverage
mkdir -p ../rmsd_qsub_coverage
mkdir -p ../padiff_qsub_coverage

for ((coverage=10; coverage<=100; coverage+=10)); do
    rm -f ../vf_qsub_coverage/*
    rm -f ../du_dn_llik_qsub_coverage/*
    rm -f ../params_qsub_coverage/old${coverage}
    rm -f ../params_qsub_coverage/new${coverage}
    rm -f ../params_qsub_coverage/best${coverage}
    rm -f ../data_coverage/* # remove data files
    rm -f ../llik_qsub_coverage/llik${coverage}
    rm -f ../rmsd_qsub_coverage/rmsd${coverage}
    rm -f ../padiff_qsub_coverage/padiff${coverage}
    
    param0=../params_qsub_coverage/old${coverage}
    cp ../../read_generation/generated/subtype2_topology0_snv/params/0.4_0.8_4.x_y_u_n $param0
    
    pa_true=../../read_generation/generated/subtype2_topology0_snv/params/0.1_0.3_3.x_y_u_n

    ./split_coverage.py ../../read_generation/generated/subtype2_topology0_snv/reads/0.1/0.3/3/coverage${coverage}_snv1000_seed1.reads 1000 100
    
    fileNum=`ls -l ../data_coverage/ | wc -l`
    fileNum=$(($fileNum - 1))

    vf0=../vf_qsub_coverage/old${coverage}
    ../bin/calc_vf_dvf_multinomial 2 0 $param0 $vf0

    qsub -N em_coverage${coverage} -sync y mstep_coverage.sh $coverage $fileNum $pa_true
done
