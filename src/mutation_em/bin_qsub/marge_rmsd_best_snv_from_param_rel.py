#!/usr/bin/env python

import sys
from subprocess import check_output
from math import sqrt

filebase = sys.argv[1]
params_filebase = sys.argv[2]
snvs = int(sys.argv[3])
maxiter = int(sys.argv[4])

bestfile = filebase + '/bests_rel'
with open(bestfile, "w") as bestf:
    for i in range(1,11,1):
        filename_i = params_filebase + '/' + str(i) + '/'
        rel_rmsd_best = 1000000.0
        rel_rmsd_best_t1 = 1000000.0
        rel_rmsd_best_t2 = 1000000.0
        rel_rmsd_best_n1 = 1000000.0
        rel_rmsd_best_n2 = 1000000.0
        rel_rmsd_best_t = 1000000.0
        rel_rmsd_best_n = 1000000.0
        rel_rmsd_best_index = 0
        for j in range(1,maxiter+1,1):
            filename_i_j = filename_i + str(j) + '/rmsd_rel'
            with open(filename_i_j, "r") as rmsd_i_j:
                l = rmsd_i_j.readline().strip().split('\t')
                rel_rmsd_t1 = float(l[0])
                rel_rmsd_t2 = float(l[1])
                rel_rmsd_n1 = float(l[2])
                rel_rmsd_n2 = float(l[3])
                rel_rmsd_t = float(l[4])
                rel_rmsd_n = float(l[5])
                rel_rmsd_t_n = float(l[6])
                
                if rel_rmsd_t_n < rel_rmsd_best:
                    rel_rmsd_best_t1 = rel_rmsd_t1
                    rel_rmsd_best_t2 = rel_rmsd_t2
                    rel_rmsd_best_n1 = rel_rmsd_n1
                    rel_rmsd_best_n2 = rel_rmsd_n2
                    rel_rmsd_best_t = rel_rmsd_t
                    rel_rmsd_best_n = rel_rmsd_n
                    rel_rmsd_best_index = j
                    rel_rmsd_best = rel_rmsd_t_n
                    
        bestf.write("%d\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n" % (snvs, rel_rmsd_best_t1, rel_rmsd_best_t2, rel_rmsd_best_n1, rel_rmsd_best_n2, rel_rmsd_best_t, rel_rmsd_best_n, rel_rmsd_best))
    bestf.close()
