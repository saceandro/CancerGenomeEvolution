#!/usr/bin/python

dic = {}
count = {}
for i in range(100):
    f = open("./mixture_cn_io/mixture_cn_accuracy_randomseed.accuracy" + str(i+1), "r")
    for line in f.readlines():
        item = line.rstrip('\n').split('\t')
        n = int(item[0])
        err = float(item[1])
        if n in dic.keys():
            dic[n] += err
            count[n] += 1
        else:
            dic[n] = err
            count[n] = 1
    f.close()

g = open("mixture_cn_accuracy_randomseed.accuracy.mean", "w")
for n in sorted(dic.keys()):
    g.write("%d\t%f\n" % (n, dic[n] / count[n]))
g.close()
