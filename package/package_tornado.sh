#!/bin/sh

# package_tornado.sh
#
# Make a tarball out of all the tornado scripts 
#

VER=`cat VERSION`
ARNAME=tornado-${VER}.tar.gz
tar -zcvf ${ARNAME} --exclude '*.pyc' --exclude '*.tar.gz' src

echo "s3cmd put --acl-public ${ARNAME} s3://tornado-emr/bin/"
