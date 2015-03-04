
#!/usr/bin/env bash

# Set location of Rail repo here
RAILHOME=~/rail

# Set output directory here -- must be on S3!
OUTPUT=s3://rail-results/geuv_ami_3_4_0b

python $RAILHOME/src prep elastic -m GEUVADIS_all_descriptive.manifest -c 20 --core-instance-bid-price 0.08 --master-instance-bid-price 0.08 --core-instance-type c3.2xlarge --master-instance-type c3.2xlarge -o $OUTPUT --ec2-key-name rail2 --do-not-check-manifest -f