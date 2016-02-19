#!/usr/bin/python

f = open("mixture_cn_accuracy_randomseed.accuracy", "r")
g = open("mixture_cn_accuracy_randomseed.accuracy.mean", "w")

dic = {}
lines = f.readlines()
num = len(lines)
for line in lines:
    item = line.split('\t').rstrip('\n')
    n = froat(item[0])
    err = froat(item[1])
    if n in dic:
        dic[n] += err
    else:
        dic[n] = 0

for n in dic.keys():
    print "%i\t%d\n" % (n, dic[n] / num)
