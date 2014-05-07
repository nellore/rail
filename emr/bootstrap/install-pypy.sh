#!/bin/bash

# install-pypy.sh
#
# Install PyPy 2.2.1; stored on s3 because bitbucket downloads are rate-limited

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

# wget --no-check-certificate -t 3 https://bitbucket.org/squeaky/portable-pypy/downloads/pypy-2.2.1-linux_x86_64-portable.tar.bz2 || { echo 'wget failed' ; exit 1; }
tar xvjf pypy-2.2.1-linux_x86_64-portable.tar.bz2 || { echo 'unzipping failed' ; exit 1; }
sudo mv pypy-2.2.1-linux_x86_64-portable /opt/pypy || { echo 'sudo mv failed' ; exit 1; }
sudo ln -s /opt/pypy/bin/pypy /usr/local/bin || { echo 'sudo ln failed' ; exit 1; }
