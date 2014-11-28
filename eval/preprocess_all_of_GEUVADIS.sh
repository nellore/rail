#!/usr/bin/env bash

# Set location of Rail repo here
RAILHOME=~/rail

# Set output directory here -- must be on S3!
OUTPUT=s3://rail-papers/ALLGEUVADISDAT

python $RAILHOME/src prep elastic -m GEUVADIS_all_descriptive.manifest -c 5 --core-instance-bid-price 0.10 --master-instance-bid-price 0.10 -o $OUTPUT --ec2-key-name rail2 --do-not-check-manifest --keep-alive --termination-protected