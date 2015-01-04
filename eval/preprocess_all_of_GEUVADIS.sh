#!/usr/bin/env bash

# Set location of Rail repo here
RAILHOME=~/rail

# Set output directory here -- must be on S3!
OUTPUT=s3://rail-west-2/prepallofgeuvadis2

python $RAILHOME/src prep elastic -m GEUVADIS_all_descriptive.manifest -c 50 --core-instance-bid-price 0.03 --master-instance-bid-price 0.03 --core-instance-type m3.large --master-instance-type m3.large -o $OUTPUT --ec2-key-name westernpacificrail --do-not-check-manifest --no-consistent-view --region us-west-2