#!/bin/zsh
./binomial_generate binomial_generate.out100 0.1 100000 100
./binomial_inference_multisites 100 binomial_generate.out100
