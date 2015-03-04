#!/bin/bash

# install-bowtie2.sh
#
# Install bowtie version 2.2.2 as of now.

set -e

export HOME=/home/hadoop
sudo wget -O/etc/yum.repos.d/s3tools.repo http://s3tools.org/repo/RHEL_6/s3tools.repo || { echo 'wget failed' ; exit 1; }
sudo yum -y install s3cmd || { echo 's3cmd installation failed' ; exit 1; }

AWS_ACCESS_ID=`grep 'fs.s3.awsAccessKeyId' ~/conf/*.xml | sed 's/.*<value>//' | sed 's/<\/value>.*//'`
AWS_ACCESS_KEY=`grep 'fs.s3.awsSecretAccessKey' ~/conf/*.xml | sed 's/.*<value>//' | sed 's/<\/value>.*//'`

cat >~/.s3cfg <<EOF
[default]
access_key = $AWS_ACCESS_ID
secret_key = $AWS_ACCESS_KEY
socket_timeout = 300
multipart_chunk_size_mb = 15
reduced_redundancy = False
send_chunk = 4096
EOF

s3cmd get s3://rail-emr/bin/bowtie2-2.2.4-linux-x86_64.zip|| { echo 's3cmd get failed' ; exit 1; }
unzip bowtie2-2.2.4-linux-x86_64.zip || { echo 'unzip failed' ; exit 1; }
sudo ln -s `pwd`/bowtie2-2.2.4/bowtie2 /usr/local/bin
sudo ln -s `pwd`/bowtie2-2.2.4/bowtie2-build /usr/local/bin