#!/bin/bash

# s3cmd_s3_tarball.sh
#
# Download a tarball from an S3 bucket and expand to given directory
#
# Arguments are:
# 1. s3:// URL to copy from
# 2. Local directory to expand archive in

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
socket_timeout = 300
multipart_chunk_size_mb = 15
reduced_redundancy = False
send_chunk = 4096
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
