#!/usr/bin/env bash

# Set location of Rail repo here
RAILHOME=~/rail

# Set output directory here -- must be on S3!
OUTPUT=s3://rail-west-2/prepallofgeuvadis

python $RAILHOME/src prep elastic -m GEUVADIS_all_descriptive.manifest -c 50 --core-instance-bid-price 0.02 --master-instance-bid-price 0.02 --core-instance-type m1.large --master-instance-type m1.large -o $OUTPUT --ec2-key-name westernpacificrail --do-not-check-manifest --no-consistent-view --region us-west-2