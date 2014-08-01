#!/usr/bin/env bash
# Set location of Rail repo here
RAILHOME=~/railclones/rail

# Set input/output bucket here -- must be on S3!
BUCKET=s3://rail-runs

python $RAILHOME/src align elastic -i s3://rail-experiments/GEUV -a hg19 -o $BUCKET/GEUVADIS_vAug.1.2014 -c 200 -m $RAILHOME/eval/GEUVADIS_all_samples.manifest