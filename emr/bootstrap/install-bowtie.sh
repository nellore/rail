#!/bin/bash

# install-bowtie.sh
#
# Install bowtie version 1.0.1 as of now.

set -e
export HOME=/home/hadoop
printf '\nexport HOME=/home/hadoop\n' >>/home/hadoop/conf/hadoop-user-env.sh

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

sudo ln -s /home/hadoop/.s3cfg /home/.s3cfg

s3cmd get s3://rail-emr/bin/bowtie-1.1.0-linux-x86_64.zip || { echo 's3cmd failed' ; exit 1; }
unzip bowtie-1.1.0-linux-x86_64.zip || { echo 'unzip failed' ; exit 1; }
sudo ln -s `pwd`/bowtie-1.1.0/bowtie /usr/local/bin
sudo ln -s `pwd`/bowtie-1.1.0/bowtie-build /usr/local/bin