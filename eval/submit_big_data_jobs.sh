#!/usr/bin/env bash

# Set location of Rail repo here
RAIL_HOME=~/railclones/rail

# Set base location of output on S3 here
OUT_BASENAME=s3://rail-experiments/geuvadis_sim

# Run Rail-RNA on 20 simulated "bioreps", pooling info across them
python $RAIL_HOME/src go elastic -a hg19 -m $RAIL_HOME/eval/geuvadis_sim.manifest -o ${OUT_BASENAME}_by_job_out -c 50 --core-instance-bid-price 0.11 --master-instance-bid-price 0.11

# Run Rail-RNA on 20 simulated "bioreps" without pooling
python $RAIL_HOME/src go elastic -a hg19 -m $RAIL_HOME/eval/geuvadis_sim.manifest -o ${OUT_BASENAME}_by_techrep_out -c 50 --core-instance-bid-price 0.11 --master-instance-bid-price 0.11 --by techrep --do-not-check-manifest