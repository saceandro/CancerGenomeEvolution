#!/usr/bin/env python

import sys

filebase = "../rmsd_qsub_snv_seed_rand/"
params_filebase = "../params_qsub_snv_seed_rand/"
bestfile = filebase + 'bests'

with open(bestfile, "w") as bestf:
    for snvs in range(1000,5001,1000):
        filebase_snvs = params_filebase + str(snvs) + '/'
        for i in range(1,11,1):
            rmsd_file = filebase_snvs + str(i) + '/1/rmsd'
            with open(rmsd_file, "r") as rmsd_f:
                l = rmsd_f.readline().strip().split('\t')
                bestf.write("%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (snvs, l[0], l[1], l[2], l[3], l[4], l[5], l[6]))
                rmsd_f.close()
    bestf.close()
