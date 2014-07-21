#!/usr/bin/env bash

# Set location of Rail repo here
RAILHOME=~/railclones/rail

# Set output directory here -- must be on S3!
OUTPUT=s3://rail-experiments/GEUV80

python $RAILHOME/src prep elastic -m GEUVADIS_80_samples.manifest -c 20 --core-instance-bid-price 0.11 --master-instance-bid-price 0.11 -o $OUTPUT --ec2-key-name rail --do-not-check-manifest