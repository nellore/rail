#!/usr/bin/env bash

# Set location of Rail repo here
RAILHOME=~/rail

# Set output directory here -- must be on S3!
OUTPUT=s3://rail-results/geuv112_ami_3_4_0e

python $RAILHOME/src prep elastic -m $RAILHOME/eval/GEUVADIS_112_staged.manifest -c 20 --max-task-attempts 10 --core-instance-type c3.2xlarge --master-instance-type c3.2xlarge --core-instance-bid-price 0.11 --master-instance-bid-price 0.11 -o $OUTPUT --ec2-key-name rail2 --do-not-check-manifest