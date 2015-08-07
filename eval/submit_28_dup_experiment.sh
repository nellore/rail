#!/usr/bin/env bash
# Use Rail-RNA v0.1.8a
# Set input/output S3 bucket here -- must be consistent with that in prep_28_dup_experiment.sh!
BUCKET=s3://rail-eu-west-1
# Get current working directory
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

rail-rna align elastic -i $BUCKET/28duppreprocessed -a hg19 -o $BUCKET/28dup_out -c 10 --master-instance-type c3.2xlarge --core-instance-type c3.2xlarge -m ${DIR}/GEUVADIS_28_with_112_sample_labels.manifest --ec2-key-name raileuw1 --region eu-west-1 --master-instance-bid-price 0.11 --core-instance-bid-price 0.11

