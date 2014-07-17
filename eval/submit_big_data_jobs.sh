#!/usr/bin/env bash

# Set location of Rail repo here
RAILHOME=~/railclones/rail

# Set base location of output on S3 here
OUTBASENAME=s3://rail-experiments/geuvadis_sim

# Run Rail-RNA on 20 simulated "bioreps", pooling info across them
python $RAILHOME/src go elastic -a hg19 -m $RAILHOME/eval/geuvadis_sim.manifest -o ${OUTBASENAME}_across_out -c 50 --core-instance-bid-price 0.11 --master-instance-bid-price 0.11 -f --ec2-key-name rail

# Run Rail-RNA on 20 simulated "bioreps" without pooling
python $RAILHOME/src go elastic -a hg19 -m $RAILHOME/eval/geuvadis_sim.manifest -o ${OUTBASENAME}_by_techrep_out -c 50 --core-instance-bid-price 0.11 --master-instance-bid-price 0.11 --by techrep --do-not-check-manifest -f --ec2-key-name rail