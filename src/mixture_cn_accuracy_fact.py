#!/usr/bin/python

from math import pow
g = open("mixture_cn_accuracy_randomseed.accuracy.mean", "w")
dic = {}
count = {}
for i in range(10):
    num = pow(2, i+1)
    f = open("./mixture_cn_io/mixture_cn_accuracy.accuracy" + str(int(num)), "r")
    for line in f.readlines():
        g.write(line)
    f.close()
g.close()
