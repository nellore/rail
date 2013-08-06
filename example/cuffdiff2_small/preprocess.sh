#!/bin/sh

if [ -z "$TORNADO_HOME" ] ; then
	echo "TORNADO_HOME is not set"
	exit 1
fi

cat cuffdiff2_small.manifest | \
python $TORNADO_HOME/src/rnawesome/preprocess.py --gzip-output
mkdir -p preprocessed_reads
mv *.gz preprocessed_reads
