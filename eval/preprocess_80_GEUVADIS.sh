#!/usr/bin/env bash

# Set location of Rail repo here
RAILHOME=~/rail

# Set output directory here -- must be on S3!
OUTPUT=s3://rail-results/GEUV80

python $RAILHOME/src prep elastic -m GEUVADIS_80_samples.manifest -c 20 --core-instance-bid-price 0.17 --master-instance-bid-price 0.17 -o $OUTPUT --ec2-key-name rail2 --do-not-check-manifest