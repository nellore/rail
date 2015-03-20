#!/usr/bin/env bash

# Set location of Rail repo here
RAILHOME=~/rail

# Set output directory here -- must be on S3!
BUCKET=s3://rail-results

python $RAILHOME/src prep elastic -m $RAILHOME/eval/GEUVADIS_112_descriptive.manifest -c 20 --core-instance-type c3.2xlarge --master-instance-type c3.2xlarge --core-instance-bid-price 0.11 --master-instance-bid-price 0.11 -o $BUCKET/112preprocessed --ec2-key-name rail2 --do-not-check-manifest
python $RAILHOME/src prep elastic -m $RAILHOME/eval/GEUVADIS_56_with_112_sample_labels.manifest -c 20 --core-instance-type c3.2xlarge --master-instance-type c3.2xlarge --core-instance-bid-price 0.11 --master-instance-bid-price 0.11 -o $BUCKET/56preprocessed --ec2-key-name rail2 --do-not-check-manifest
python $RAILHOME/src prep elastic -m $RAILHOME/eval/GEUVADIS_28_with_112_sample_labels.manifest -c 20 --core-instance-type c3.2xlarge --master-instance-type c3.2xlarge --core-instance-bid-price 0.11 --master-instance-bid-price 0.11 -o $BUCKET/28preprocessed --ec2-key-name rail2 --do-not-check-manifest