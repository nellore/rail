#!/bin/sh

python ../../src/simulate/simple.py \
	--output-prefix eg1 \
	--fasta ../fasta/virus.fa \
	--seed 77 \
	--num-nucs 100000
