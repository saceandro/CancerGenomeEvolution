#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N em_qsub_snv
#$ -e ../log_qsub_snv/em.err
#$ -o ../log_qsub_snv/em.log

u_lower=$1
em_max_iter=$2
grad_desc_max_iter=$3

source ~/.zshrc
rm -f ../log_qsub_snv/esteplog
rm -f ../log_qsub_snv/esteperr

mkdir -p ../log_qsub_snv
mkdir -p ../llik_qsub_snv
mkdir -p ../vf_qsub_snv
mkdir -p ../params_qsub_snv
mkdir -p ../rmsd_qsub_snv
mkdir -p ../padiff_qsub_snv

#snvs=100
#snvs=19694
# for ((snvs=10; snvs<=100; snvs*=10)); do
#     mkdir -p ../data_snv/${snvs}
#     mkdir -p ../du_dn_llik_qsub_snv/${snvs}

#     rm -f ../vf_qsub_snv/old${snvs}
#     rm -f ../vf_qsub_snv/new${snvs}
#     rm -f ../du_dn_llik_qsub_snv/${snvs}/*
#     rm -f ../params_qsub_snv/old${snvs}
#     rm -f ../params_qsub_snv/new${snvs}
#     rm -f ../params_qsub_snv/best${snvs}
#     rm -f ../data_snv/${snvs}/* # remove data files
#     rm -f ../llik_qsub_snv/llik${snvs}
#     rm -f ../rmsd_qsub_snv/rmsd${snvs}
#     rm -f ../padiff_qsub_snv/padiff${snvs}
#     rm -f ../log_qsub_snv/log${snvs}
#     rm -f ../log_qsub_snv/err${snvs}
    
#     param0=../params_qsub_snv/old${snvs}
# #    cp ../../read_generation/generated/subtype2_topology0_snv/params/0.4_0.8_4.x_y_u_n $param0
#     cp ../../read_generation/generated/subtype2_topology0_snv/params/0.1_0.3_3.x_y_u_n $param0
     
#     pa_true=../../read_generation/generated/subtype2_topology0_snv/params/0.1_0.3_3.x_y_u_n

# #    ./split_snv.py ../../read_generation/generated/subtype2_topology0/reads/0.1/0.3/3/coverage100_bp30000000_seed1.reads $snvs 100
#     ./split_snv.py ../../read_generation/generated/subtype2_topology0_snv/reads/0.1/0.3/3/coverage100_snv${snvs}_seed1.reads $snvs 100
    
#     fileNum=`ls -l ../data_snv/${snvs} | wc -l`
#     fileNum=$(($fileNum - 1))

#     vf0=../vf_qsub_snv/old${snvs}
#     ../bin/calc_vf_dvf_multinomial 2 0 $param0 $vf0

#     qsub -N em_snv${snvs} mstep_snv.sh $snvs $fileNum $pa_true $u_lower $em_max_iter $grad_desc_max_iter
# done

for ((snvs=1000; snvs<=10000; snvs*=10)); do
    mkdir -p ../data_snv/${snvs}
    mkdir -p ../du_dn_llik_qsub_snv/${snvs}

    rm -f ../vf_qsub_snv/old${snvs}
    rm -f ../vf_qsub_snv/new${snvs}
    rm -f ../du_dn_llik_qsub_snv/${snvs}/*
    rm -f ../params_qsub_snv/old${snvs}
    rm -f ../params_qsub_snv/new${snvs}
    rm -f ../params_qsub_snv/best${snvs}
    rm -f ../data_snv/${snvs}/* # remove data files
    rm -f ../llik_qsub_snv/llik${snvs}
    rm -f ../rmsd_qsub_snv/rmsd${snvs}
    rm -f ../padiff_qsub_snv/padiff${snvs}
    rm -f ../log_qsub_snv/log${snvs}
    rm -f ../log_qsub_snv/err${snvs}
    
    param0=../params_qsub_snv/old${snvs}
#    cp ../../read_generation/generated/subtype2_topology0_snv/params/0.4_0.8_4.x_y_u_n $param0
    cp ../../read_generation/generated/subtype2_topology0_snv/params/0.1_0.3_3.x_y_u_n $param0
     
    pa_true=../../read_generation/generated/subtype2_topology0_snv/params/0.1_0.3_3.x_y_u_n

#    ./split_snv.py ../../read_generation/generated/subtype2_topology0/reads/0.1/0.3/3/coverage100_bp30000000_seed1.reads $snvs 100
    ./split_snv.py ../../read_generation/generated/subtype2_topology0_snv/reads/0.1/0.3/3/coverage100_snv${snvs}_seed1.reads $snvs 100
    
    fileNum=`ls -l ../data_snv/${snvs} | wc -l`
    fileNum=$(($fileNum - 1))

    vf0=../vf_qsub_snv/old${snvs}
    ../bin/calc_vf_dvf_multinomial 2 0 $param0 $vf0

    qsub -N em_snv${snvs} -sync y mstep_snv.sh $snvs $fileNum $pa_true $u_lower $em_max_iter $grad_desc_max_iter
done
