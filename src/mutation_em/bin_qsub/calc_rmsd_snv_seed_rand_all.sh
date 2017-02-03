#!/bin/env zsh

for ((snvs=1000; snvs<=5000; snvs+=1000)); do
    ./calc_rmsd_snv_seed_rand_param.py ../params_qsub_snv_seed_rand/$snvs ../../read_generation/generated/subtype2_topology0_snv/params/0.5_0.4_3.t_n 1
done
