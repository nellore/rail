#!/bin/sh

# package_tornado.sh
#
# Make a tarball out of all the tornado scripts 
#

VER=`cat VERSION`
ARNAME=tornado_${VER}.tar.gz
cd src
tar -zcvf ${ARNAME} --exclude '*.pyc' --exclude '*.tar.gz' .
mv ${ARNAME} ..
cd ..
