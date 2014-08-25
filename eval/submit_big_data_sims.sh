#!/usr/bin/env bash

# Set location of Rail repo here
RAILHOME=~/railclones/rail

# Set base location of output on S3 here
OUTBASENAME=s3://rail-experiments/geuvadis_sim

# Run Rail-RNA on 20 simulated "bioreps", pooling info across them
python $RAILHOME/src align elastic -a hg19 -m $RAILHOME/eval/geuvadis_sim.manifest -i s3://rail-experiments/GEUVSIM20.intermediate/preprocess/push -o s3://rail-experiments/GEUVADISSIM20REPS2 -c 30 --ec2-key-name rail --core-instance-bid-price 0.20 --master-instance-bid-price 0.20 --master-instance-type c3.2xlarge --core-instance-type c3.2xlarge -f