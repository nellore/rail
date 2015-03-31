
#!/usr/bin/env bash

# Set location of Rail repo here
RAILHOME=~/rail

# Set output directory here -- must be on S3!
OUTPUT=s3://rail-west-2/geuvprepped

python $RAILHOME/src prep elastic -m GEUVADIS_all_descriptive_staged.manifest -c 84 --core-instance-bid-price 0.11 --master-instance-bid-price 0.11 --core-instance-type c3.2xlarge --master-instance-type c3.2xlarge -o $OUTPUT --ec2-key-name westernpacificrail -f --do-not-check-manifest