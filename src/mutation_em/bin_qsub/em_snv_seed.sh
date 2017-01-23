#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N em_qsub_snv_seed
#$ -e ../log_qsub_snv_seed/em.err
#$ -o ../log_qsub_snv_seed/em.log

u_lower=$1
em_max_iter=$2
grad_desc_max_iter=$3

source ~/.zshrc
rm -f ../log_qsub_snv_seed/esteplog
rm -f ../log_qsub_snv_seed/esteperr

for ((snvs=1000; snvs<=9000; snvs+=2000)); do
    echo "snv: " $snvs
    for ((seed=1; seed<=10; ++seed)); do
        echo "seed: " $seed
        mkdir -p ../log_qsub_snv_seed/${snvs}
        mkdir -p ../llik_qsub_snv_seed/${snvs}
        mkdir -p ../vf_qsub_snv_seed/${snvs}/${seed}
        mkdir -p ../params_qsub_snv_seed/${snvs}/${seed}
        mkdir -p ../rmsd_qsub_snv_seed/${snvs}
        mkdir -p ../padiff_qsub_snv_seed/${snvs}
        mkdir -p ../data_snv_seed/${snvs}/${seed}
        mkdir -p ../du_dn_llik_qsub_snv_seed/${snvs}/${seed}

        rm -f ../vf_qsub_snv_seed/${snvs}/${seed}/old
        rm -f ../vf_qsub_snv_seed/${snvs}/${seed}/new
        rm -f ../du_dn_llik_qsub_snv_seed/${snvs}/${seed}/*
        rm -f ../params_qsub_snv_seed/${snvs}/${seed}/old
        rm -f ../params_qsub_snv_seed/${snvs}/${seed}/new
        rm -f ../params_qsub_snv_seed/${snvs}/${seed}/best
        rm -f ../params_qsub_snv_seed/${snvs}/${seed}/log
        rm -f ../data_snv_seed/${snvs}/${seed}/* # remove data files
        rm -f ../llik_qsub_snv_seed/${snvs}/${seed}
        rm -f ../rmsd_qsub_snv_seed/${snvs}/${seed}
        rm -f ../padiff_qsub_snv_seed/${snvs}/${seed}
        rm -f ../log_qsub_snv_seed/${snvs}/${seed}.log
        rm -f ../log_qsub_snv_seed/${snvs}/${seed}.err
        

        u1=0.$(($RANDOM % 9 + 1))
        u2=0.$(($RANDOM % 9 + 1))
        n=$(($RANDOM % 9 + 1))
        u1u2n=${u1}_${u2}_${n}
        echo $u1u2n

        param0=../params_qsub_snv_seed/${snvs}/${seed}/old
        cp ../../read_generation/generated/subtype2_topology0_snv/params/${u1u2n}.x_y_u_n $param0
        
        pa_true=../../read_generation/generated/subtype2_topology0_snv/params/0.5_0.4_3.x_y_u_n

        ./split_snv_seed.py ../../read_generation/generated/subtype2_topology0_snv/reads/0.5/0.4/3/coverage100_snv${snvs}_seed${seed}.reads $snvs $seed 100
        
        fileNum=`ls -l ../data_snv_seed/${snvs}/${seed} | wc -l`
        fileNum=$(($fileNum - 1))

        vf0=../vf_qsub_snv_seed/${snvs}/${seed}/old
        ../bin/calc_vf_dvf_multinomial 2 0 $param0 $vf0

        qsub -N em_snv${snvs}_seed${seed} -sync y mstep_snv_seed.sh $snvs $seed $fileNum $pa_true $u_lower $em_max_iter $grad_desc_max_iter
    done
    echo "-------------------------------------------------------------------------"
done
