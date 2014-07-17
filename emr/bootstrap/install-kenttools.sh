#!/bin/bash

# install-kenttools.sh
#
# Install some of the more useful tools from Jim Kent for manipulating bed,
# wig, bigBed and bigWig files.

set -e

wget -O- -q http://s3tools.org/repo/deb-all/stable/s3tools.key | sudo apt-key add -
sudo wget -O/etc/apt/sources.list.d/s3tools.list http://s3tools.org/repo/deb-all/stable/s3tools.list
sudo apt-get update && sudo apt-get install s3cmd

AWS_ACCESS_ID=`grep 'fs.s3.awsAccessKeyId' $HOME/conf/*.xml | sed 's/.*<value>//' | sed 's/<\/value>.*//'`
AWS_ACCESS_KEY=`grep 'fs.s3.awsSecretAccessKey' $HOME/conf/*.xml | sed 's/.*<value>//' | sed 's/<\/value>.*//'`

cat >$HOME/.s3cfg <<EOF
[default]
access_key = $AWS_ACCESS_ID
secret_key = $AWS_ACCESS_KEY
EOF

# Change to destination directory
mkdir -p $1
cd $1

if [ `uname -m` != "x86_64" ] ; then
	echo "Not a 64-bit platform!"
	exit 1
fi

s3cmd get s3://rail-emr/bin/bedToBigBed || { echo 's3cmd get failed' ; exit 1; }
chmod a+x bedToBigBed
