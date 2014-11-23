#!/usr/bin/env bash
# Set location of Rail repo here
RAILHOME=~/rail

# Set input/output bucket here -- must be on S3!
BUCKET=s3://rail-paper

# Submit 28 at 10, 20, and 40 cores
python $RAILHOME/src align elastic -i $BUCKET/GEUV28 -a hg19 -o $BUCKET/GEUVADIS_28_datasets_10_c3.2xlarge -c 10 --master-instance-type c3.2xlarge --core-instance-type c3.2xlarge -m $RAILHOME/eval/GEUVADIS_28_with_112_sample_labels.manifest --ec2-key-name rail2
python $RAILHOME/src align elastic -i $BUCKET/GEUV28 -a hg19 -o $BUCKET/GEUVADIS_28_datasets_20_c3.2xlarge -c 20 --master-instance-type c3.2xlarge --core-instance-type c3.2xlarge -m $RAILHOME/eval/GEUVADIS_28_with_112_sample_labels.manifest --ec2-key-name rail2
#python $RAILHOME/src align elastic -i $BUCKET/GEUV28 -a hg19 -o $BUCKET/GEUVADIS_28_datasets_40_c3.2xlarge -c 40 --master-instance-type c3.2xlarge --core-instance-type c3.2xlarge -m $RAILHOME/eval/GEUVADIS_28_with_112_sample_labels.manifest --ec2-key-name rail2

# Submit 56 and 112 at 40 cores
python $RAILHOME/src align elastic -i $BUCKET/GEUV56 -a hg19 -o $BUCKET/GEUVADIS_56_datasets_40_c3.2xlarge -c 40 --master-instance-type c3.2xlarge --core-instance-type c3.2xlarge -m $RAILHOME/eval/GEUVADIS_56_with_112_sample_labels.manifest --ec2-key-name rail2
python $RAILHOME/src align elastic -i $BUCKET/GEUV112 -a hg19 -o $BUCKET/GEUVADIS_112_datasets_40_c3.2xlarge -c 40 --master-instance-type c3.2xlarge --core-instance-type c3.2xlarge -m $RAILHOME/eval/GEUVADIS_112.manifest --ec2-key-name rail2