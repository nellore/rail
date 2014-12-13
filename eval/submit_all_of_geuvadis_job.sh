#!/usr/bin/env bash
# Set location of Rail repo here
RAILHOME=~/rail

# Set input/output bucket here -- must be on S3!
BUCKET=s3://rail-papers

python $RAILHOME/src align elastic -i s3://rail-papers/geuvadisprepped -a hg19 -o $BUCKET/GEUVADIS_v12.11.2014b -c 42 -m $RAILHOME/eval/GEUVADIS_all_descriptive.manifest --master-instance-type c3.8xlarge --core-instance-type c3.8xlarge --master-instance-bid-price 0.27 --core-instance-bid-price 0.27 --ec2-key-name rail2