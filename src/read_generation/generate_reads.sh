#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N generate_reads_poisson

./generate_reads_poisson 4 1000000 1000 2 30000000000 1 generated/4_2_original.u_n generated/4_2_original_poisson_1000000_1000_1_cell1000000.reads 1> generated/4_2_original_poisson_1000000_1000_1_cell1000000.log 2> generated/4_2_original_poisson_1000000_1000_1_cell1000000.err
