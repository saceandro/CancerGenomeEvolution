#!/usr/bin/env python

import sys
from subprocess import check_output
from math import sqrt

filebase = "../rmsd_qsub_snv_seed_rand_param/"
filebase_u1 = "../rmsd_qsub_snv_seed_rand_param/0.5/"
u2_1 = sys.argv[1]
u2_2 = sys.argv[2]
# u2_3 = sys.argv[3]
marged_bestfile = filebase + 'bests'
bestfile1 = filebase_u1 + u2_1 + '/3/3000/bests'
bestfile2 = filebase_u1 + u2_2 + '/3/3000/bests'
# bestfile3 = filebase_u1 + u2_3 + '/3/3000/bests'

with open(marged_bestfile, "w") as bestf:
    with open(bestfile1, "r") as bestf1:
        for i in range(10):
            l = bestf1.readline().strip().split('\t')
            bestf.write("%s\t%.2f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (u2_1, float(u2_1) * 0.5, l[1], l[2], l[3], l[4], l[5], l[6], l[7]))
        bestf1.close()
    with open(bestfile2, "r") as bestf2:
        for i in range(10):
            l = bestf2.readline().strip().split('\t')
            bestf.write("%s\t%.2f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (u2_2, float(u2_2) * 0.5, l[1], l[2], l[3], l[4], l[5], l[6], l[7]))
        bestf2.close()
    # with open(bestfile3, "r") as bestf3:
    #     for i in range(10):
    #         l = bestf3.readline().strip().split('\t')
    #         bestf.write("%s\t%.2f\t%s\t%s\t%s\n" % (u2_3, float(u2_3) * 0.5, l[1], l[2], l[3]))
    #     bestf3.close()
    bestf.close()
