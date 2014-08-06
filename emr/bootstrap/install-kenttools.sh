#!/bin/bash

# install-kenttools.sh
#
# Install some of the more useful tools from Jim Kent for manipulating bed,
# wig, bigBed and bigWig files.

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

# Change to destination directory
mkdir -p $1
cd $1

if [ `uname -m` != "x86_64" ] ; then
	echo "Not a 64-bit platform!"
	exit 1
fi

s3cmd get s3://rail-emr/bin/bedGraphToBigWig || { echo 's3cmd get failed' ; exit 1; }
chmod a+x bedGraphToBigWig
sudo ln -s /home/hadoop/bedGraphToBigWig /bin/bedGraphToBigWig
