#!/usr/bin/env bash

# Set location of Rail repo here
RAILHOME=~/railclones/rail

# Set output directory here -- must be on S3!
OUTPUT=s3://rail-results/GEUV100

python $RAILHOME/src prep elastic -m GEUVADIS_100_samples.manifest -c 20 --core-instance-bid-price 0.17 --master-instance-bid-price 0.17 -o $OUTPUT --ec2-key-name rail --do-not-check-manifest