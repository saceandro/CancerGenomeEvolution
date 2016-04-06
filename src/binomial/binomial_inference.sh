#!/bin/sh
./binomial_generate binomial_generate.out 0.1 1000 1000
./binomial_inference binomial_generate.out 1000 binomial_inference.out
