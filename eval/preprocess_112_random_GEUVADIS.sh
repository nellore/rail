#!/usr/bin/env bash

# Set location of Rail repo here
RAILHOME=~/rail

# Set output directory here -- must be on S3!
OUTPUT=s3://rail-paper/GEUV100

python $RAILHOME/src prep elastic -m $RAILHOME/eval/GEUVADIS_112.manifest -c 20 --core-instance-type m3.large --master-instance-type m3.large --core-instance-bid-price 0.10 --master-instance-bid-price 0.10 -o $OUTPUT --ec2-key-name rail2