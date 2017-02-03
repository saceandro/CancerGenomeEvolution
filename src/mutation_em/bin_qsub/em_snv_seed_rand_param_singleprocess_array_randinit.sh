#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N em_qsub_snv_seed_rand_param
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

snvs=3000

u1=0.5
u2=0.4
n=1

pa_true=../../read_generation/generated/subtype2_topology0_snv/params/${u1}_${u2}_${n}.x_y_u_n

for ((seed=1; seed<=10; ++seed)); do
    echo $seed
    mkdir -p ../log_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}
    mkdir -p ../llik_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}
    mkdir -p ../rmsd_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}
    mkdir -p ../padiff_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}
    mkdir -p ../data_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}
    rm -f ../data_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/* # remove data files

    read_data=../../read_generation/generated/subtype2_topology0_snv/reads/${u1}/${u2}/${n}/coverage100_snv${snvs}_seed${seed}.reads
    ./split_snv_seed_rand_param.py $read_data $snvs $seed 1 $u1 $u2 $n # 1 split

    for ((iter=1; iter<=10; ++iter)); do
        echo $iter
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

        ./calc_reverse_param_rand $subtypes $topology $param0 $iter
#        cp ../../read_generation/generated/subtype2_topology0_snv/params/${initu1u2n}.x_y_u_n $param0

        vf0=../vf_qsub_snv_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${seed}/${iter}/old
        ../bin/calc_vf_dvf_multinomial 2 0 $param0 $vf0
        echo "-----------------------------------------------------------------------------------------------"
    done
    echo "================================================================================================="
done

fileNum=1
qsub -N em_snv${snvs}_seed${seed}_rand_param${u1}_${u2}_${n} -t 1-100:1 mstep_snv_seed_rand_param_singleprocess_array.sh $u1 $u2 $n $snvs $fileNum $pa_true $u_lower $em_max_iter $grad_desc_max_iter
