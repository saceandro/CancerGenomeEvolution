#!/usr/bin/python

f = open("mixture_cn_accuracy_randomseed.accuracy", "r")
g = open("mixture_cn_accuracy_randomseed.accuracy.mean", "w")

dic = {}
count = {}
lines = f.readlines()
num = len(lines)
for line in lines:
    item = line.rstrip('\n').split('\t')
    n = int(item[0])
    err = float(item[1])
    if n in dic.keys():
        dic[n] += err
        count[n] += 1
    else:
        dic[n] = err
        count[n] = 1

for n in sorted(dic.keys()):
    g.write("%d\t%f\n" % (n, dic[n] / count[n]))
