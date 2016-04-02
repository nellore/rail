#!/usr/bin/env bash

# Set output bucket here
BUCKET=s3://rail-eu-west-1
# Current dir
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

rail-rna prep elastic -m ${DIR}/GEUVADIS_28_dup_experiment.manifest -c 20 --core-instance-type c3.2xlarge --master-instance-type c3.2xlarge --core-instance-bid-price 0.11 --master-instance-bid-price 0.11 -o $BUCKET/28dupprep --ec2-key-name raileuw1 --do-not-check-manifest