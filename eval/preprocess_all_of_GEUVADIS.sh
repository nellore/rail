#!/usr/bin/env bash

# Set location of Rail repo here
RAILHOME=~/railclones/rail

# Set output directory here -- must be on S3!
OUTPUT=s3://rail-experiments/PREPROCESSEDGEUVADIS

python $RAILHOME/src prep elastic -m GEUVADIS_all_samples.manifest -c 20 --core-instance-bid-price 0.15 --master-instance-bid-price 0.15 -o $OUTPUT --ec2-key-name rail --do-not-check-manifest