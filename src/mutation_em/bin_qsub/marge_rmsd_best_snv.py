#!/usr/bin/env python

import sys
from subprocess import check_output
from math import sqrt

filebase = sys.argv[1]
snvs = int(sys.argv[2])
maxiter = int(sys.argv[3])

bestfile = filebase + '/bests'
with open(bestfile, "w") as bestf:
    for i in range(1,11,1):
        filename_i = filebase + '/' + str(i) + '/'
        rmsd_best = 1000000.0
        rmsd_best_u = 1000000.0
        rmsd_best_t = 1000000.0
        rmsd_best_n = 1000000.0
        rmsd_best_index = 0
        for j in range(1,maxiter+1,1):
            filename_i_j = filename_i + str(j)
            l = check_output(['tail', '-n 1', filename_i_j]).strip().split('\t')
            rmsd_u = float(l[1])
            rmsd_t = float(l[2])
            rmsd_n = float(l[3])
            rmsd_t_n = sqrt((rmsd_t*rmsd_t + rmsd_n*rmsd_n)/2.0)
            if rmsd_t_n < rmsd_best:
                rmsd_best = rmsd_t_n
                rmsd_best_u = rmsd_u
                rmsd_best_t = rmsd_t
                rmsd_best_n = rmsd_n
                rmsd_best_index = j
        bestf.write("%d\t%.10e\t%.10e\t%.10e\n" % (snvs, rmsd_best_u, rmsd_best_t, rmsd_best_n))
    bestf.close()
