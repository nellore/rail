#!/bin/bash

# s3cmd_s3.sh
#
# Download a file from an S3 bucket to given directory.  Optionally rename it.
#
# Arguments are:
# 1. s3:// URL to copy from
# 2. Local directory to copy to
# 3. If specified, name to rename file to

set -e

sudo wget -O/etc/yum.repos.d/s3tools.repo http://s3tools.org/repo/RHEL_6/s3tools.repo || { echo 'wget failed' ; exit 1; }
sudo yum install s3cmd || { echo 's3cmd installation failed' ; exit 1; }

AWS_ACCESS_ID=`grep 'fs.s3.awsAccessKeyId' $HOME/conf/*.xml | sed 's/.*<value>//' | sed 's/<\/value>.*//'`
AWS_ACCESS_KEY=`grep 'fs.s3.awsSecretAccessKey' $HOME/conf/*.xml | sed 's/.*<value>//' | sed 's/<\/value>.*//'`

cat >$HOME/.s3cfg <<EOF
[default]
access_key = $AWS_ACCESS_ID
secret_key = $AWS_ACCESS_KEY
socket_timeout = 300
multipart_chunk_size_mb = 15
reduced_redundancy = False
send_chunk = 4096
EOF

mkdir -p ${2}
cd ${2}
fn=`basename ${1}`
s3cmd get ${1} . || { echo 's3cmd get failed' ; exit 1; }
if [ -n "${3}" ] ; then
	mv $fn ${3} || true
fi
