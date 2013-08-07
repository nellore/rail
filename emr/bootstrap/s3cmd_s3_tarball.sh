#!/bin/bash

# s3cmd_s3_tarball.sh
#
# Download a tarball from an S3 bucket and expand to given directory
#
# Arguments are:
# 1. s3:// URL to copy from
# 2. Local directory to expand archive in

set -e

sudo apt-get --yes install s3cmd || { echo 'apt-get failed' ; exit 1; }

AWS_ACCESS_ID=`grep 'fs.s3.awsAccessKeyId' $HOME/conf/*.xml | sed 's/.*<value>//' | sed 's/<\/value>.*//'`
AWS_ACCESS_KEY=`grep 'fs.s3.awsSecretAccessKey' $HOME/conf/*.xml | sed 's/.*<value>//' | sed 's/<\/value>.*//'`

cat >$HOME/.s3cfg <<EOF
[default]
access_key = $AWS_ACCESS_ID
secret_key = $AWS_ACCESS_KEY
EOF

mkdir -p $1
cd $1
shift

for i in $* ; do
	s3cmd get $i - | gzip -dc | tar xvf - &
	PIDS="$PIDS $!"
done

wait $PIDS
if [ $? -ne 0 ] ; then
	echo 'one or more s3cmds failed'
	exit 1
fi
