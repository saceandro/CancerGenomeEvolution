#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N em_qsub_snv_seed_rand_param
#$ -e ../log_qsub_snv_seed_rand_param/em.err
#$ -o ../log_qsub_snv_seed_rand_param/em.log

source ~/.zshrc

u_lower=$1
em_max_iter=$2
grad_desc_max_iter=$3
snvs=$4
seed=$5
u1=$6
u2=$7
n=$8

pa_true=../../read_generation/generated/subtype2_topology0_snv/params/${u1}_${u2}_${n}.x_y_u_n
read_data=../../read_generation/generated/subtype2_topology0_snv/reads/${u1}/${u2}/${n}/coverage100_snv${snvs}_seed${seed}.reads

echo "snv: " $snvs

    echo "seed: " $seed
    mkdir -p ../log_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}
    mkdir -p ../llik_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}
    mkdir -p ../rmsd_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}
    mkdir -p ../padiff_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}
    mkdir -p ../data_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}
#    rm -f ../data_snv_seed_rand_param/${snvs}/${seed}/* # remove data files

    ./split_snv_seed_rand_param.py $read_data $snvs $seed 100 $u1 $u2 $n
    fileNum=`ls -l ../data_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed} | wc -l`
    fileNum=$(($fileNum - 1))

    for ((iter=1; iter<=10; ++iter)); do
        echo "iter: " $iter
        mkdir -p ../vf_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}
        mkdir -p ../params_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}
        mkdir -p ../du_dn_llik_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}
        
        rm -f ../vf_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}/old
        rm -f ../vf_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}/new
        rm -f ../du_dn_llik_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}/*
        rm -f ../params_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}/old
        rm -f ../params_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}/new
        rm -f ../params_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}/best
        rm -f ../params_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}/log
        rm -f ../llik_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}
        rm -f ../rmsd_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}
        rm -f ../padiff_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}
        rm -f ../log_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}.log
        rm -f ../log_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}.err

        initu1=0.$(($RANDOM % 9 + 1))
        initu2=0.$(($RANDOM % 9 + 1))
        initn=$(($RANDOM % 9 + 1))
        initu1u2n=${initu1}_${initu2}_${initn}
        echo $initu1u2n

        param0=../params_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}/old
        cp ../../read_generation/generated/subtype2_topology0_snv/params/${initu1u2n}.x_y_u_n $param0

        vf0=../vf_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}/old
        ../bin/calc_vf_dvf_multinomial 2 0 $param0 $vf0

        qsub -N em_snv${snvs}_seed${seed}_rand_param_iter${iter} -sync y mstep_snv_seed_rand_param.sh $u1 $u2 $n $snvs $seed $iter $fileNum $pa_true $u_lower $em_max_iter $grad_desc_max_iter
    done
    echo "-------------------------------------------------------------------------"
