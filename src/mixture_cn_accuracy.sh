#!/bin/sh
./mixture_cn_generate mixture_cn_generate.out100 100000 100 4 kappa.txt
./mixture_cn_accuracy 4 4 100 mixture_cn_generate.out100 mixture_cn_accuracy.out100 kappa.txt mixture_cn_accuracy.accuracy100
