#!/usr/bin/env bash

# Set input/output bucket here -- must be on S3!
BUCKET=s3://rail-eu-west-1
# Current dir
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# Chance --ec2-key-name or remove it below
rail-rna align elastic -i $BUCKET/geuv667prepped -a hg19 -o $BUCKET/GEUVADIS_v07.14.2015 -c 60 -m ${DIR}/GEUVADIS_all_descriptive_revised.manifest --master-instance-type c3.8xlarge --core-instance-type c3.8xlarge --master-instance-bid-price 0.40 --core-instance-bid-price 0.40 --ec2-key-name raileuw1 -f --max-task-attempts 10 --region eu-west-1