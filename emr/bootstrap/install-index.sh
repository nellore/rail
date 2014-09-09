#!/usr/bin/env bash

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

mkdir -p $1
cd $1

fn=`basename $2`

cat >.download_idx.py <<EOF
"""
.download_idx.py

Checks if file size hasn't changed while downloading idx, and if it hasn't,
kills subprocess and restarts it.
"""

import sys
import time
import os
import subprocess

tries = 0
url = sys.argv[1]
filename = sys.argv[2]
while tries < 5:
    break_outer_loop = False
    s3cmd_process \
        = subprocess.Popen(['s3cmd', 'get', url, './', '-f'])
    time.sleep(1)
    last_check_time = time.time()
    try:
        last_size = os.path.getsize(filename)
    except OSError:
        last_size = 0
    while s3cmd_process.poll() is None:
        now_time = time.time()
        if now_time - last_check_time > 120:
            try:
                new_size = os.path.getsize(filename)
            except OSError:
                new_size = 0
            if new_size == last_size:
                # Download stalled
                break_outer_loop = True
                break
            else:
                last_size = new_size
            last_check_time = now_time
            time.sleep(1)
    if break_outer_loop:
        tries += 1
        s3cmd_process.kill()
        try:
            os.remove(filename)
        except OSError:
            pass
        time.sleep(2)
        continue
    if s3cmd_process.poll() > 0:
        if tries > 5:
            raise RuntimeError(('Non-zero exitlevel %d from s3cmd '
                                'get command')  % (
                                                s3cmd_process.poll()
                                            ))
        else:
            try:
                os.remove(filename)
            except OSError:
                pass
            tries += 1
            time.sleep(2)
            continue                       
    break
if tries > 5:
    raise RuntimeError('Could not download file from S3 '
                       'within 5 tries.')

EOF

python .download_idx.py $2 $fn
tar xvzf $fn
cd index
ln -s genome.4.bt2 genome.4.ebwt
cd ..