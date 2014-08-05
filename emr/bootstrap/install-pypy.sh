#!/bin/bash

# install-pypy.sh
#
# Install PyPy 2.2.1; stored on s3 because bitbucket downloads are rate-limited

set -e

sudo wget -O/etc/yum.repos.d/s3tools.repo http://s3tools.org/repo/RHEL_6/s3tools.repo || { echo 'wget failed' ; exit 1; }
sudo yum install s3cmd || { echo 's3cmd installation failed' ; exit 1; }

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

s3cmd get s3://rail-emr/bin/pypy-2.2.1-linux_x86_64-portable.tar.bz2 || { echo 's3cmd failed' ; exit 1; }

# wget --no-check-certificate -t 3 https://bitbucket.org/squeaky/portable-pypy/downloads/pypy-2.2.1-linux_x86_64-portable.tar.bz2 || { echo 'wget failed' ; exit 1; }
tar xvjf pypy-2.2.1-linux_x86_64-portable.tar.bz2 || { echo 'unzipping failed' ; exit 1; }
sudo mv pypy-2.2.1-linux_x86_64-portable /opt/pypy || { echo 'sudo mv failed' ; exit 1; }
sudo ln -s /opt/pypy/bin/pypy /usr/local/bin || { echo 'sudo ln failed' ; exit 1; }
