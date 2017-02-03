#!/usr/bin/env python

import sys

filebase = "../rmsd_qsub_snv_seed_rand/"
params_filebase = "../params_qsub_snv_seed_rand/"
bestfile = filebase + 'bests_rel'

with open(bestfile, "w") as bestf:
    for snvs in range(1000,5001,1000):
        filebase_snvs = params_filebase + str(snvs) + '/'
        for i in range(1,11,1):
            rmsd_rel_file = filebase_snvs + str(i) + '/1/rmsd_rel'
            with open(rmsd_rel_file, "r") as rmsd_rel_f:
                l = rmsd_rel_f.readline().strip().split('\t')
                bestf.write("%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (snvs, l[0], l[1], l[2], l[3], l[4], l[5], l[6]))
                rmsd_rel_f.close()
    bestf.close()
