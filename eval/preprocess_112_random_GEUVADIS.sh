#!/usr/bin/env bash

# Set location of Rail repo here
RAILHOME=~/rail

# Set output directory here -- must be on S3!
OUTPUT=s3://rail-paper/GEUV112

python $RAILHOME/src prep elastic -m $RAILHOME/eval/GEUVADIS_112.manifest -c 20 --core-instance-type m3.xlarge --master-instance-type m3.xlarge --core-instance-bid-price 0.15 --master-instance-bid-price 0.15 -o $OUTPUT --ec2-key-name rail2 --do-not-check-manifest