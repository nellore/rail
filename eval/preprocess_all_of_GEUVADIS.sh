
#!/usr/bin/env bash

# Set location of Rail repo here
RAILHOME=~/rail

# Set output directory here -- must be on S3!
OUTPUT=s3://rail-eu-west-1/geuvprepped

python $RAILHOME/src prep elastic -m GEUVADIS_all_descriptive_staged.manifest -c 1 -o $OUTPUT -f --region eu-west-1 --do-not-check-manifest --ec2-key-name raileuw1