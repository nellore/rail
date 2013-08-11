#!/bin/sh

# package_tornado.sh
#
# Make a tarball out of all the tornado scripts 
#

VER=`cat VERSION`
ARNAME=tornado-${VER}.tar.gz
cd src && make clean && cd ..

make -C hadoop

# Make lib directory for the multiplefiles.jar file
mkdir -p lib
cp hadoop/*.jar lib

tar -zcvf ${ARNAME} --exclude '*.pyc' --exclude '*.tar.gz' src lib

echo "s3cmd put --acl-public ${ARNAME} s3://tornado-emr/bin/"
