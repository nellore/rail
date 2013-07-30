#!/bin/bash

# s3cmd_s3_tarball.sh
#
# Download a tarball from an S3 bucket and expand to given directory
#
# Arguments are:
# 1. s3:// URL to copy from
# 2. Local directory to expand archive in

set -e

sudo apt-get --yes install s3cmd

AWS_ACCESS_ID=`grep 'fs.s3.awsAccessKeyId' $HOME/conf/*.xml | sed 's/.*<value>//' | sed 's/<\/value>.*//'`
AWS_ACCESS_KEY=`grep 'fs.s3.awsSecretAccessKey' $HOME/conf/*.xml | sed 's/.*<value>//' | sed 's/<\/value>.*//'`

cat >$HOME/.s3cfg <<EOF
[default]
access_key = $AWS_ACCESS_ID
secret_key = $AWS_ACCESS_KEY
EOF

fn=`basename $1`
s3cmd get ${1} .
mkdir -p ${2}
tar -xzf $fn -C ${2}
rm -f $fn
