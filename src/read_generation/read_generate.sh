#!/bin/env zsh
#$ -S /usr/local/bin/zsh
#$ -cwd
#$ -N generate_reads_allinherited

# ./generate_reads_allinherited 4 10000 1000 2 1 generated/4_2_original.u_n generated/4_2_original_10000_1000_1_cell100000.reads 1> generated/4_2_original_10000_1000_1_cell100000.log 2> generated/4_2_original_10000_1000_1_cell100000.err
# ./generate_reads_allinherited 4 10000 1000 2 1 generated/4_2_original.u_n generated/4_2_original_10000_1000_1_cell100000_2.reads 1> generated/4_2_original_10000_1000_1_cell100000_2.log 2> generated/4_2_original_10000_1000_1_cell100000_2.err
# ./generate_reads_allinherited 4 1000000 1000 2 1 generated/4_2_original.u_n generated/4_2_original_1000000_1000_1_cell100000.reads 1> generated/4_2_original_1000000_1000_1_cell100000.log 2> generated/4_2_original_1000000_1000_1_cell100000.err
# ./generate_reads_allinherited 4 1000000 1000 2 1 generated/4_2_original.u_n generated/4_2_original_1000000_1000_1_cell100000000.reads 1> generated/4_2_original_1000000_1000_1_cell100000000.log 2> generated/4_2_original_1000000_1000_1_cell100000000.err
./generate_reads_allinherited 4 1000000 1000 2 1 generated/4_2_original.u_n generated/4_2_original_1000000_1000_1_cell1000000.reads 1> generated/4_2_original_1000000_1000_1_cell1000000.log 2> generated/4_2_original_1000000_1000_1_cell1000000.err
