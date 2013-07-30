#!/bin/bash

# s3cmd_s3_tarball.sh
#
# Download a tarball from an S3 bucket and expand to given directory
#
# Arguments are:
# 1. Bucket to copy from
# 2. Subdir to copy from
# 3. Local directory to expand archive in

set -e

sudo apt-get --yes install s3cmd

AWS_ACCESS_ID=`grep 'fs.s3.awsAccessKeyId' $HOME/conf/*.xml | sed 's/.*<value>//' | sed 's/<\/value>.*//'`
AWS_ACCESS_KEY=`grep 'fs.s3.awsSecretAccessKey' $HOME/conf/*.xml | sed 's/.*<value>//' | sed 's/<\/value>.*//'`

cat >$HOME/.s3cfg <<EOF
[default]
access_key = $AWS_ACCESS_ID
secret_key = $AWS_ACCESS_KEY
EOF

fn=`basename $2`
s3cmd get s3://${1}/${2} .
mkdir -p ${3}
tar -xzf $fn -C ${3}
rm -f $fn
