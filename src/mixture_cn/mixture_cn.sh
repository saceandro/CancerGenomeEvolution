#!/bin/sh
./mixture_cn_generate mixture_cn_generate.out10000 10000 10000 4 kappa.txt
./mixture_cn 4 4 10000 mixture_cn_generate.out10000 mixture_cn.out10000
