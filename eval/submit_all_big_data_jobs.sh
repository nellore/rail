#!/usr/bin/env bash
# Set location of Rail repo here
RAILHOME=~/railclones/rail

# Set input/output bucket here -- must be on S3!
BUCKET=s3://rail-experiments

# Submit 40 at 25, 50, and 100 cores
python $RAILHOME/src align elastic -i $BUCKET/GEUV40 -a hg19 -o $BUCKET/GEUVADIS40_25cores_out_2 -c 25 -m $RAILHOME/eval/GEUVADIS_40_samples.manifest
python $RAILHOME/src align elastic -i $BUCKET/GEUV40 -a hg19 -o $BUCKET/GEUVADIS40_50cores_out_2 -c 50 -m $RAILHOME/eval/GEUVADIS_40_samples.manifest
python $RAILHOME/src align elastic -i $BUCKET/GEUV40 -a hg19 -o $BUCKET/GEUVADIS40_100cores_out_2 -c 100 -m $RAILHOME/eval/GEUVADIS_40_samples.manifest

# Submit 20 and 80 at 50 cores
#python $RAILHOME/src align elastic -i $BUCKET/GEUV20 -a hg19 -o $BUCKET/GEUVADIS20_50cores_out -c 50 -m $RAILHOME/eval/GEUVADIS_20_samples.manifest
python $RAILHOME/src align elastic -i $BUCKET/GEUV80 -a hg19 -o $BUCKET/GEUVADIS80_50cores_out_2 -c 50 -m $RAILHOME/eval/GEUVADIS_80_samples.manifest