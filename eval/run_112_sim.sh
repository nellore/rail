#!/usr/bin/env bash
## A script that:
## (1) Compresses the FASTQs and stages them on S3
## (2) Creates a manifest file for the FASTQs
## (3) Runs Rail-RNA on the 112 datasets on Elastic MapReduce
## First command-line parameter: directory with fastqs
## Second command-line parameter: where on S3 to stage fastq.gzs
## Third command-line parameter: path to manifest file to write
## Requires s3cmd
FASTQDIR=$1
S3STAGED=$2
S3DEST=$3
MANIFEST=$4
cd $FASTQDIR
ls *.fastq | python -c "import sys; for line in sys.stdin:
    sample_name = list(line.strip().rpartition('_')[0].partition('-'))
    sample_name[0] = sample_name[0] + '_sim'
    print '\t'.join([${S3STAGED}/line.strip(), '0', ''.join(sample_name)])"
#for i in *.fastq; do cat $i | gzip | s3cmd put - $S3DEST/$i.gz; done