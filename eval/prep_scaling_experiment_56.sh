#!/usr/bin/env bash

# Set output directory here -- must be on S3!
BUCKET=s3://rail-eu-west-1
# Current dir
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

rail-rna prep elastic -m ${DIR}/GEUVADIS_56_with_112_sample_labels.manifest -c 20 --core-instance-type c3.2xlarge --master-instance-type c3.2xlarge --core-instance-bid-price 0.13 --master-instance-bid-price 0.13 -o $BUCKET/56preprocessed --ec2-key-name raileuw1 --do-not-check-manifest