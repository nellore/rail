#!/usr/bin/env bash
# Set location of Rail repo here
RAILHOME=~/railclones/rail

# Set input/output bucket here -- must be on S3!
BUCKET=s3://rail-experiments

python $RAILHOME/src align elastic -i $BUCKET/GEUV -a hg19 -o $BUCKET/GEUV_out -c 300 -m $RAILHOME/eval/GEUVADIS_all_samples.manifest