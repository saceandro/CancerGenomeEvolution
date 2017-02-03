#!/usr/bin/env python

import sys

f_name = sys.argv[1]
file_len = int(sys.argv[2])
split_num = int(sys.argv[3])
lines = 0
remainder = file_len % split_num
if remainder==0:
    lines = file_len / split_num
else:
    lines = file_len / split_num + 1

rest = file_len

if lines > 101:
    with open(f_name, "r") as f:
        for i in range(split_num):
            if rest > 0:
                with open("../data/" + str(i), "w") as g:
                    writelines = min(lines,rest)
                    for j in range(writelines):
                        g.write(f.readline())
                    g.close()
                rest -= writelines
            else:
                break
        f.close()

else:
    with open(f_name, "r") as f:
        for i in range(101):
            if rest > 0:
                with open("../data/" + str(i), "w") as g:
                    writelines = min(100,rest)
                    for j in range(writelines):
                        g.write(f.readline())
                    g.close()
                rest -= writelines
            else:
                break
        f.close()


    # for i in range(split_num-1):
    #     with open("data/" + str(i), "w") as g:
    #         for j in range(lines):
    #             g.write(f.readline())
    #         g.close()
    # with open("data/" + str(split_num-1), "w") as g:
    #     for j in range(lines+remainder):
    #         g.write(f.readline())
    #     g.close()
    # f.close()
