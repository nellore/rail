#!/usr/bin/env bash
# Set location of Rail repo here
RAILHOME=~/rail

# Set input/output bucket here -- must be on S3 and consistent with that in prep_scaling_experiments.sh!
BUCKET=s3://rail-results

# Submit 28 at 10, 20, and 40 cores
python $RAILHOME/src align elastic -i $BUCKET/28preprocessed -a hg19 -o $BUCKET/28out.10cores -c 10 --master-instance-type c3.2xlarge --core-instance-type c3.2xlarge -m $RAILHOME/eval/GEUVADIS_28_with_112_sample_labels.manifest --ec2-key-name rail2 --master-instance-bid-price 0.11 --core-instance-bid-price 0.11
python $RAILHOME/src align elastic -i $BUCKET/28preprocessed -a hg19 -o $BUCKET/28out.20cores -c 20 --master-instance-type c3.2xlarge --core-instance-type c3.2xlarge -m $RAILHOME/eval/GEUVADIS_28_with_112_sample_labels.manifest --ec2-key-name rail2 --master-instance-bid-price 0.11 --core-instance-bid-price 0.11
python $RAILHOME/src align elastic -i $BUCKET/28preprocessed -a hg19 -o $BUCKET/28out.40cores -c 40 --master-instance-type c3.2xlarge --core-instance-type c3.2xlarge -m $RAILHOME/eval/GEUVADIS_28_with_112_sample_labels.manifest --ec2-key-name rail2 --master-instance-bid-price 0.11 --core-instance-bid-price 0.11

# Submit 56 and 112 at 40 cores
python $RAILHOME/src align elastic -i $BUCKET/56preprocessed -a hg19 -o $BUCKET/56out.40cores -c 40 --master-instance-type c3.2xlarge --core-instance-type c3.2xlarge -m $RAILHOME/eval/GEUVADIS_56_with_112_sample_labels.manifest --ec2-key-name rail2 --master-instance-bid-price 0.11 --core-instance-bid-price 0.11
python $RAILHOME/src align elastic -i $BUCKET/128preprocessed -a hg19 -o $BUCKET/112out.40cores -c 40 --master-instance-type c3.2xlarge --core-instance-type c3.2xlarge -m $RAILHOME/eval/GEUVADIS_112.manifest --ec2-key-name rail2 --master-instance-bid-price 0.11 --core-instance-bid-price 0.11