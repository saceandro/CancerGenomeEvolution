#!/usr/bin/env python

import sys
from math import sqrt

filebase = sys.argv[1]
truefile = sys.argv[2]
maxiter = int(sys.argv[3])

true_t1=0
true_t2=0
true_n1=0
true_n2=0
with open(truefile, "r") as truef:
    line = truef.readline().strip().split('\t')
    true_t1 = float(line[1])
    true_t2 = float(line[2])
    lline = truef.readline().strip().split('\t')
    true_n1 = float(lline[1])
    true_n2 = float(lline[2])
    truef.close()
    
for i in range(1,11,1):
    filename_i = filebase + '/' + str(i) + '/'
    for j in range(1,maxiter+1,1):
        filename_i_j = filename_i + str(j) + '/best'
        rmsdfile = filename_i + str(j) + '/rmsd'
        rel_rmsdfile = filename_i + str(j) + '/rmsd_rel'

        with open(filename_i_j, "r") as bestparam_f:
            bestparam_f.readline()
            bestparam_f.readline()
            l = bestparam_f.readline().strip().split('\t')
            t1 = float(l[0])
            t2 = t1 * float(l[1])
            ll = bestparam_f.readline().strip().split('\t')
            n1 = float(ll[1])
            n2 = float(ll[2])
            
            abs_t1 = abs(t1 - true_t1)
            abs_t2 = abs(t2 - true_t2)
            abs_n1 = abs(n1 - true_n1)
            abs_n2 = abs(n2 - true_n2)
            abs_t = sqrt((abs_t1 ** 2 + abs_t2 ** 2) / 2.0)
            abs_n = sqrt((abs_n1 ** 2 + abs_n2 ** 2) / 2.0)
            abs_t_n = sqrt((abs_t ** 2 + abs_n ** 2) / 2.0)
            
            rel_t1 = abs_t1 / true_t1
            rel_t2 = abs_t2 / true_t2
            rel_n1 = abs_n1 / true_n1
            rel_n2 = abs_n2 / true_n2
            rel_t = sqrt((rel_t1 ** 2 + rel_t2 ** 2) / 2.0)
            rel_n = sqrt((rel_n1 ** 2 + rel_n2 ** 2) / 2.0)
            rel_t_n = sqrt((rel_t ** 2 + rel_n ** 2) / 2.0)
            
            with open(rmsdfile, "w") as rmsdf:
#                rmsdf.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("abs_t1", "abs_t2", "abs_n1", "abs_n2", "abs_t", "abs_n", "abs_t_n"))
                rmsdf.write("%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n" % (abs_t1, abs_t2, abs_n1, abs_n2, abs_t, abs_n, abs_t_n))
                rmsdf.close()
                    
            with open(rel_rmsdfile, "w") as rel_rmsdf:
#                rel_rmsdf.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("rel_t1", "rel_t2", "rel_n1", "rel_n2", "rel_t", "rel_n", "rel_t_n"))
                rel_rmsdf.write("%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n" % (rel_t1, rel_t2, rel_n1, rel_n2, rel_t, rel_n, rel_t_n))
                rel_rmsdf.close()
                
            bestparam_f.close()
