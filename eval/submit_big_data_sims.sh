#!/usr/bin/env bash

# Set location of Rail repo here
RAILHOME=~/railclones/rail

# Set base location of output on S3 here
OUTBASENAME=s3://rail-experiments/geuvadis_sim

# Run Rail-RNA on 20 simulated "bioreps", pooling info across them
python $RAILHOME/src align elastic -a hg19 -m $RAILHOME/eval/geuvadis_sim.manifest -i s3://rail-experiments/testhadoop2.1.intermediate/preprocess/push -o s3://rail-experiments/testhadoop2onbioreps3 -c 50 --master-instance-bid-price 0.13 --core-instance-bid-price 0.13 --ec2-key-name rail --master-instance-type c3.2xlarge --core-instance-type c3.2xlarge --keep-alive --termination-protected -f