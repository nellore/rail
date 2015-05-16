#!/bin/bash

# get_s3_tarball.sh
#
# Used for bootstrap actions where a publicly-visible tarball is being
# downloaded from an S3 bucket and expanded to a given directory
#
# Arguments are:
# 1. Bucket to copy from
# 2. Subdir to copy from
# 3. Local directory to expand archive in

set -e
export HOME=/home/hadoop
fn=`basename $2`
wget -S -T 10 -t 5 http://${1}.s3.amazonaws.com/${2} || { echo 'wget failed' ; exit 1; }
mkdir -p ${3}
tar -xzf $fn -C ${3} || { echo 'untar failed' ; exit 1; }
rm -f $fn
