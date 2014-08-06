#!/usr/bin/env bash

set -e

export HOME=/home/hadoop

AWS_ACCESS_ID=`grep 'fs.s3.awsAccessKeyId' $HOME/conf/*.xml | sed 's/.*<value>//' | sed 's/<\/value>.*//'`
AWS_ACCESS_KEY=`grep 'fs.s3.awsSecretAccessKey' $HOME/conf/*.xml | sed 's/.*<value>//' | sed 's/<\/value>.*//'`

wget https://s3.amazonaws.com/aws-cli/awscli-bundle.zip || { echo 'wget failed' ; exit 1; }
unzip awscli-bundle.zip
sudo ./awscli-bundle/install -i /usr/local/aws -b /usr/local/bin/aws

mkdir -p $HOME/.aws
cat >~/.aws/config <<EOF
[default]
aws_access_key_id = $AWS_ACCESS_ID
aws_secret_access_key = $AWS_ACCESS_KEY
EOF

mkdir -p $1
cd $1

fn=`basename $2`

aws s3 cp $2 ./
tar xzvf $fn