#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N em_qsub_coverage_seed_rand_param
#$ -e ../log_qsub_coverage_seed_rand_param/em.err
#$ -o ../log_qsub_coverage_seed_rand_param/em.log
#$ -l s_vmem=1G,mem_req=1G

source ~/.zshrc

u_lower=$1
em_max_iter=$2
grad_desc_max_iter=$3

mkdir -p ../log_qsub_coverage_seed_rand_param
rm -f ../log_qsub_coverage_seed_rand_param/em.log
rm -f ../log_qsub_coverage_seed_rand_param/em.err
rm -f ../log_qsub_coverage_seed_rand_param/esteplog
rm -f ../log_qsub_coverage_seed_rand_param/esteperr

snvs=3000

u1=0.5
u2=0.4
n=3

pa_true=../../read_generation/generated/subtype2_topology0_snv/params/${u1}_${u2}_${n}.x_y_u_n

for ((coverage=10; coverage<=100; coverage*=10)); do
    for ((seed=1; seed<=10; ++seed)); do
        echo $seed
        mkdir -p ../log_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}
        mkdir -p ../llik_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}
        mkdir -p ../rmsd_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}
        mkdir -p ../padiff_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}
        mkdir -p ../data_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}
        rm -f ../data_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/* # remove data files

        read_data=../../read_generation/generated/subtype2_topology0_snv/reads/${u1}/${u2}/${n}/coverage${coverage}_snv${snvs}_seed${seed}.reads
        cp $read_data ../data_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/0

        for ((iter=1; iter<=5; ++iter)); do
            echo $iter
            mkdir -p ../vf_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}
            mkdir -p ../params_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}
            mkdir -p ../du_dn_llik_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}

            rm -f ../vf_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}/old
            rm -f ../vf_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}/new
            rm -f ../du_dn_llik_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}/*
            rm -f ../params_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}/old
            rm -f ../params_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}/new
            rm -f ../params_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}/best
            rm -f ../params_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}/log
            rm -f ../llik_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}
            rm -f ../rmsd_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}
            rm -f ../padiff_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}
            rm -f ../log_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}.log
            rm -f ../log_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}.err

            param0=../params_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}/old

            ./calc_reverse_param_rand_n 2 0 $param0 $iter 0.4 0
            vf0=../vf_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}/old
            ../bin/calc_vf_dvf_multinomial 2 0 $param0 $vf0
            echo "-----------------------------------------------------------------------------------------------"
        done
        echo "================================================================================================="
    done
done

for ((coverage=50; coverage<=200; coverage*=4)); do
    for ((seed=1; seed<=10; ++seed)); do
        echo $seed
        mkdir -p ../log_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}
        mkdir -p ../llik_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}
        mkdir -p ../rmsd_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}
        mkdir -p ../padiff_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}
        mkdir -p ../data_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}
        rm -f ../data_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/* # remove data files

        read_data=../../read_generation/generated/subtype2_topology0_snv/reads/${u1}/${u2}/${n}/coverage${coverage}_snv${snvs}_seed${seed}.reads
        cp $read_data ../data_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/0

        for ((iter=1; iter<=5; ++iter)); do
            echo $iter
            mkdir -p ../vf_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}
            mkdir -p ../params_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}
            mkdir -p ../du_dn_llik_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}

            rm -f ../vf_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}/old
            rm -f ../vf_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}/new
            rm -f ../du_dn_llik_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}/*
            rm -f ../params_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}/old
            rm -f ../params_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}/new
            rm -f ../params_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}/best
            rm -f ../params_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}/log
            rm -f ../llik_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}
            rm -f ../rmsd_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}
            rm -f ../padiff_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}
            rm -f ../log_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}.log
            rm -f ../log_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}.err

            param0=../params_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}/old

            ./calc_reverse_param_rand_n 2 0 $param0 $iter 0.4 0
            vf0=../vf_qsub_coverage_seed_rand_param/${u1}/${u2}/${n}/${snvs}/${coverage}/${seed}/${iter}/old
            ../bin/calc_vf_dvf_multinomial 2 0 $param0 $vf0
            echo "-----------------------------------------------------------------------------------------------"
        done
        echo "================================================================================================="
    done
done

fileNum=1
qsub -N em_coverage_seed_rand_param${u1}_${u2}_${n} -t 1-200:1 mstep_coverage_seed_rand_param_singleprocess_array.sh $u1 $u2 $n $snvs $fileNum $pa_true $u_lower $em_max_iter $grad_desc_max_iter
