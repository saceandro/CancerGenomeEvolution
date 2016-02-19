#!/bin/zsh

qsub -t 1-100:1 mixture_cn_accuracy_array_randomseed.sh
./mixture_cn_accuracy_array_randomseed.py
