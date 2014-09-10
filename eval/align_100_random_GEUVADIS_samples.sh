#!/usr/bin/env bash
# Set location of Rail repo here
RAILHOME=~/railclones/rail

# Set input/output bucket here -- must be on S3!
BUCKET=s3://rail-experiments

python $RAILHOME/src align elastic -i $BUCKET/GEUV100 -a hg19 -o $BUCKET/GEUV100random_out -c 40 --master-instance-type c3.2xlarge --core-instance-type c3.2xlarge -m $RAILHOME/eval/GEUVADIS_100_samples.manifest --ec2-key-name rail