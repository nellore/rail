
#!/usr/bin/env bash

# Set location of Rail repo here
RAILHOME=~/rail

# Set output directory here -- must be on S3!
OUTPUT=s3://rail-eu-west-1/geuvprepped

python $RAILHOME/src prep elastic -m GEUVADIS_all_descriptive_staged.manifest -c 2 --core-instance-bid-price 0.35 --master-instance-bid-price 0.35 --core-instance-type c3.2xlarge --master-instance-type c3.2xlarge -o $OUTPUT --no-consistent-view -f --region eu-west-1 --do-not-check-manifest