#!/usr/bin/env bash
# Set location of Rail repo here
RAILHOME=~/railclones/rail

# Set input/output bucket here -- must be on S3!
BUCKET=s3://rail-experiments

python $RAILHOME/src align elastic -i $BUCKET/GEUV -a hg19 -o $BUCKET/GEUVADIS_v08.01.2014 -c 150 -m $RAILHOME/eval/GEUVADIS_all_samples.manifest --ec2-key-name rail