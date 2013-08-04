#!/bin/bash

# install-tornado.sh
#
# Install SRA toolkit, including fastq-dump.  Version 2.1.7 as of now.

set -e

AWS_ACCESS_ID=`grep 'fs.s3.awsAccessKeyId' $HOME/conf/*.xml | sed 's/.*<value>//' | sed 's/<\/value>.*//'`
AWS_ACCESS_KEY=`grep 'fs.s3.awsSecretAccessKey' $HOME/conf/*.xml | sed 's/.*<value>//' | sed 's/<\/value>.*//'`

cat >$HOME/.s3cfg <<EOF
[default]
access_key = $AWS_ACCESS_ID
secret_key = $AWS_ACCESS_KEY
EOF

fn=`basename $1`
s3cmd get ${1} . || { echo 's3cmd failed' ; exit 1; }

mkdir -p ${2}
tar -xzf $fn -C ${2} || { echo 'untar failed' ; exit 1; }
rm -f $fn
make -C ${2}/src || { echo 'make failed' ; exit 1; }
