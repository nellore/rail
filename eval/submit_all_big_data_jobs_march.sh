#!/usr/bin/env bash
# Set location of Rail repo here
RAILHOME=~/rail

# Set input/output bucket here -- must be on S3!
BUCKET=s3://rail-results

# Submit 28 at 10, 20, and 40 cores
#python $RAILHOME/src align elastic -i $BUCKET/GEUV28 -a hg19 -o $BUCKET/GEUVADIS_28_datasets_10_c3.2xlarge.revised -c 10 --master-instance-type c3.2xlarge --core-instance-type c3.2xlarge -m $RAILHOME/eval/GEUVADIS_28_with_112_sample_labels.manifest --ec2-key-name rail2
#python $RAILHOME/src align elastic -i $BUCKET/GEUV28 -a hg19 -o $BUCKET/GEUVADIS_28_datasets_20_c3.2xlarge.revised.2 -c 20 --master-instance-type c3.2xlarge --core-instance-type c3.2xlarge -m $RAILHOME/eval/GEUVADIS_28_with_112_sample_labels.manifest --ec2-key-name rail2
#python $RAILHOME/src align elastic -i $BUCKET/GEUV28 -a hg19 -o $BUCKET/GEUVADIS_28_datasets_40_c3.2xlarge.revised -c 40 --master-instance-type c3.2xlarge --core-instance-type c3.2xlarge -m $RAILHOME/eval/GEUVADIS_28_with_112_sample_labels.manifest --ec2-key-name rail2

# Submit 56 and 112 at 20 cores
#python $RAILHOME/src align elastic -i $BUCKET/GEUV56 -a hg19 -o $BUCKET/GEUVADIS_56_datasets_20_c3.2xlarge.revised.2 -c 20 --master-instance-type c3.2xlarge --core-instance-type c3.2xlarge -m $RAILHOME/eval/GEUVADIS_56_with_112_sample_labels.manifest --ec2-key-name rail2
python $RAILHOME/src align elastic -i $BUCKET/geuv112_ami_3_4_0e -a hg19 -o $BUCKET/geuv112_ami_3_4_0_out2 -c 50 --master-instance-type c3.2xlarge --core-instance-type c3.2xlarge -m $RAILHOME/eval/GEUVADIS_112.manifest --ec2-key-name rail2 --max-task-attempts 10 --master-instance-bid-price 0.11 --core-instance-bid-price 0.11