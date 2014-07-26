#!/usr/bin/env bash

# Set location of Rail repo here
RAILHOME=~/railclones/rail

# Set base location of output on S3 here
OUTBASENAME=s3://rail-experiments/geuvadis_sim

# Run Rail-RNA on 20 simulated "bioreps", pooling info across them
python $RAILHOME/src go elastic -a hg19 -m $RAILHOME/eval/geuvadis_sim.manifest -o ${OUTBASENAME}_allbioreps_pres -c 40 --core-instance-bid-price 0.11 --master-instance-bid-price 0.11 --master-instance-type c3.2xlarge --core-instance-type c3.2xlarge -f