#!/bin/sh

# package_rail.sh
#
# Make a tarball out of all the Rail scripts 
#

VER=`cat VERSION`
ARNAME=rail-rna-${VER}.tar.gz

make -C hadoop

# Make lib directory for the hacked hadoop-streaming jar
rm -rf lib
mkdir -p lib
cp hadoop/out/*.jar lib

tar -zcvf ${ARNAME} --exclude '*.pyc' --exclude '*.tar.gz' src lib

echo "s3cmd put --acl-public ${ARNAME} s3://rail-emr/bin/"
