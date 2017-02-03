#!/usr/bin/env python

import sys

f_name = sys.argv[1]
file_len = int(sys.argv[2])
seed = int(sys.argv[3])
split_num = int(sys.argv[4])
u1 = sys.argv[5]
u2 = sys.argv[6]
n = sys.argv[7]
lines = 0
remainder = file_len % split_num
if remainder==0:
    lines = file_len / split_num
else:
    lines = file_len / split_num + 1

rest = file_len

if lines > 501:
    with open(f_name, "r") as f:
        for i in range(split_num):
            if rest > 0:
                with open("../data_snv_seed_rand_param/" + u1 + "/" + u2 + "/" + n + "/" + str(file_len) + "/" + str(seed) + "/" + str(i), "w") as g:
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
        for i in range(501):
            if rest > 0:
                with open("../data_snv_seed_rand_param/" + u1 + "/" + u2 + "/" + n + "/" + str(file_len) + "/" + str(seed) + "/" + str(i), "w") as g:
                    writelines = min(500,rest)
                    for j in range(writelines):
                        g.write(f.readline())
                    g.close()
                rest -= writelines
            else:
                break
        f.close()


    # for i in range(split_num-1):
    #     with open("data_snv_seed_rand_param/" + str(i), "w") as g:
    #         for j in range(lines):
    #             g.write(f.readline())
    #         g.close()
    # with open("data_snv_seed_rand_param/" + str(split_num-1), "w") as g:
    #     for j in range(lines+remainder):
    #         g.write(f.readline())
    #     g.close()
    # f.close()
