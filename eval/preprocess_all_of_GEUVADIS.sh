#!/usr/bin/env bash

# Set location of Rail repo here
RAILHOME=~/railclones/rail

# Set output directory here -- must be on S3!
OUTPUT=s3://rail-experiments/ALLOFGEUVADIS

python $RAILHOME/src prep elastic -m GEUVADIS_all_samples.manifest -c 50 -o $OUTPUT