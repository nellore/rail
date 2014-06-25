#!/bin/bash

# install-kenttools.sh
#
# Install some of the more useful tools from Jim Kent for manipulating bed,
# wig, bigBed and bigWig files.

set -e

sudo apt-get --yes install s3cmd || { echo 'apt-get failed' ; exit 1; }

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
