
#!/usr/bin/env bash

# Set location of Rail repo here
RAILHOME=~/rail

# Set output directory here -- must be on S3!
OUTPUT=s3://rail-results/geuv_ami_3_4_0d

python $RAILHOME/src prep elastic -m GEUVADIS_all_descriptive_staged.manifest -c 84 --core-instance-bid-price 0.11 --master-instance-bid-price 0.11 --core-instance-type c3.2xlarge --master-instance-type c3.2xlarge -o $OUTPUT --ec2-key-name rail2 -f --do-not-check-manifest