#!/usr/bin/env bash

# package_rail.sh
#
# Make a tarball out of all the Rail scripts 
#

cd $(dirname "${BASH_SOURCE[0]}")
VER=`cat ../src/rna/utils/version.py | grep version_number | cut -f 2 -d "'"`
ARNAME=rail-rna-${VER}.tar.gz

# Nothing to make now; just copy jars
# make -C hadoop

# Make lib directory for the hacked hadoop-streaming jar
# rm -rf lib
# mkdir -p lib
cd ..
tar czvf ${ARNAME} --exclude '*.pyc' --exclude '*.tar.gz' src lib

echo "This should be executed from Rail's root directory!"
echo "s3cmd put --acl-public ${ARNAME} s3://rail-emr/bin/"
