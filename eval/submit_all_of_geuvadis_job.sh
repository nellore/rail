#!/usr/bin/env bash
# Set location of Rail repo here
RAILHOME=~/railclones/rail

# Set input/output bucket here -- must be on S3!
BUCKET=s3://rail-paper

python $RAILHOME/src align elastic -i $BUCKET/PREPROCESSEDGEUVADIS -a hg19 -o $BUCKET/GEUVADIS_v09.05.2014 -c 50 -m $RAILHOME/eval/GEUVADIS_all_samples.manifest --master-instance-type c3.8xlarge --core-instance-type c3.8xlarge --search-filter mild --ec2-key-name rail2