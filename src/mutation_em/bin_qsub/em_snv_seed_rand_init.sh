#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N em_qsub_snv_seed_rand
#$ -e ../log_qsub_snv_seed_rand/em.err
#$ -o ../log_qsub_snv_seed_rand/em.log

source ~/.zshrc

u_lower=$1
em_max_iter=$2
grad_desc_max_iter=$3
snvs=$4
iter=$5

u1=$6
u2=$7
n=$8

pa_true=../../read_generation/generated/subtype2_topology0_snv/params/0.5_0.4_3.x_y_u_n

echo "snv: " $snvs
for ((seed=1; seed<=10; ++seed)); do
    echo "seed: " $seed
    mkdir -p ../log_qsub_snv_seed_rand/${snvs}/${seed}
    mkdir -p ../llik_qsub_snv_seed_rand/${snvs}/${seed}
    mkdir -p ../rmsd_qsub_snv_seed_rand/${snvs}/${seed}
    mkdir -p ../padiff_qsub_snv_seed_rand/${snvs}/${seed}
    mkdir -p ../data_snv_seed_rand/${snvs}/${seed}
#    rm -f ../data_snv_seed_rand/${snvs}/${seed}/* # remove data files

    ./split_snv_seed_rand.py ../../read_generation/generated/subtype2_topology0_snv/reads/0.5/0.4/3/coverage100_snv${snvs}_seed${seed}.reads $snvs $seed 100
    fileNum=`ls -l ../data_snv_seed_rand/${snvs}/${seed} | wc -l`
    fileNum=$(($fileNum - 1))

#    for ((iter=1; iter<=5; ++iter)); do
        echo "iter: " $iter
        mkdir -p ../vf_qsub_snv_seed_rand/${snvs}/${seed}/${iter}
        mkdir -p ../params_qsub_snv_seed_rand/${snvs}/${seed}/${iter}
        mkdir -p ../du_dn_llik_qsub_snv_seed_rand/${snvs}/${seed}/${iter}
        
        rm -f ../vf_qsub_snv_seed_rand/${snvs}/${seed}/${iter}/old
        rm -f ../vf_qsub_snv_seed_rand/${snvs}/${seed}/${iter}/new
        rm -f ../du_dn_llik_qsub_snv_seed_rand/${snvs}/${seed}/${iter}/*
        rm -f ../params_qsub_snv_seed_rand/${snvs}/${seed}/${iter}/old
        rm -f ../params_qsub_snv_seed_rand/${snvs}/${seed}/${iter}/new
        rm -f ../params_qsub_snv_seed_rand/${snvs}/${seed}/${iter}/best
        rm -f ../params_qsub_snv_seed_rand/${snvs}/${seed}/${iter}/log
        rm -f ../llik_qsub_snv_seed_rand/${snvs}/${seed}/${iter}
        rm -f ../rmsd_qsub_snv_seed_rand/${snvs}/${seed}/${iter}
        rm -f ../padiff_qsub_snv_seed_rand/${snvs}/${seed}/${iter}
        rm -f ../log_qsub_snv_seed_rand/${snvs}/${seed}/${iter}.log
        rm -f ../log_qsub_snv_seed_rand/${snvs}/${seed}/${iter}.err

        u1u2n=${u1}_${u2}_${n}
        echo $u1u2n

        param0=../params_qsub_snv_seed_rand/${snvs}/${seed}/${iter}/old
        cp ../../read_generation/generated/subtype2_topology0_snv/params/${u1u2n}.x_y_u_n $param0

        vf0=../vf_qsub_snv_seed_rand/${snvs}/${seed}/${iter}/old
        ../bin/calc_vf_dvf_multinomial 2 0 $param0 $vf0

        qsub -N em_snv${snvs}_seed${seed}_iter${iter} -sync y mstep_snv_seed_rand.sh $snvs $seed $iter $fileNum $pa_true $u_lower $em_max_iter $grad_desc_max_iter
#    done
#    echo "-------------------------------------------------------------------------"
done
echo "========================================================================="
