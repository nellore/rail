#!/bin/sh

if [ -z "$RAIL_HOME" ] ; then
	echo "RAIL_HOME is not set"
	exit 1
fi

cat cuffdiff2_small.manifest | \
python $RAIL_HOME/src/rail-rna/preprocess.py --gzip-output
mkdir -p preprocessed_reads
mv *.gz preprocessed_reads
