#!/usr/bin/env bash
# Set location of Rail repo here
RAILHOME=~/rail

# Set input/output bucket here -- must be on S3!
BUCKET=s3://rail-results

python $RAILHOME/src align elastic -i s3://rail-results/geuv_ami_3_4_0d -a hg19 -o $BUCKET/GEUVADIS_v3.27.2015a -c 60 -m $RAILHOME/eval/GEUVADIS_all_descriptive_staged.manifest --master-instance-type c3.8xlarge --core-instance-type c3.8xlarge --master-instance-bid-price 0.30 --core-instance-bid-price 0.30 --ec2-key-name rail2 --max-task-attempts 10 --json