#!/bin/sh

cat cuffdiff2_small.manifest | \
python $TORNADO_HOME/src/rnawesome/preprocess.py
mkdir -p preprocessed_reads
mv *.gz preprocessed_reads
