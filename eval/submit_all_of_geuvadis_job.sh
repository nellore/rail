#!/usr/bin/env bash
# Set location of Rail repo here
RAILHOME=~/rail

# Set input/output bucket here -- must be on S3!
BUCKET=s3://rail-eu-west-1

python $RAILHOME/src align elastic -i s3://rail-eu-west-1/geuvprepped -a hg19 -o $BUCKET/GEUVADIS_v4.04.2015 -c 60 -m $RAILHOME/eval/GEUVADIS_all_descriptive_staged.manifest --master-instance-type c3.8xlarge --core-instance-type c3.8xlarge --master-instance-bid-price 0.35 --core-instance-bid-price 0.35 --ec2-key-name raileuwest1 --no-consistent-view -f --max-task-attempts 10 --region eu-west-1