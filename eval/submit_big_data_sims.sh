#!/usr/bin/env bash

# Set location of Rail repo here
RAILHOME=~/railclones/rail

# Set base location of output on S3 here
OUTBASENAME=s3://rail-experiments/geuvadis_sim

# Run Rail-RNA on 20 simulated "bioreps", pooling info across them
python $RAILHOME/src go elastic -a hg19 -m $RAILHOME/eval/geuvadis_sim.manifest -o s3://rail-experiments/testhadoop2withhdfs --intermediate hdfs:///intermediates -c 20 --ec2-key-name rail --core-instance-bid-price 0.13 --master-instance-bid-price 0.13 --master-instance-type c3.2xlarge --core-instance-type c3.2xlarge -f