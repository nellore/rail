#!/usr/bin/env bash

# package_rail.sh
#
# Make a tarball out of all the Rail scripts 
#

ORIGINAL=$(pwd)
PACKDIR=$(dirname "${BASH_SOURCE[0]}")
cd $PACKDIR
cd ..
VER=`cat src/version.py | grep version_number | cut -f 2 -d "'"`
ARNAME=rail-rna-${VER}.tar.gz
rm -rf $ARNAME
tar czvf ${ARNAME} --exclude '*.pyc' --exclude '*.tar.gz' --exclude '.DS_Store' src lib
mv $ARNAME $PACKDIR/$ARNAME
echo "Run this to upload to S3."
echo "s3cmd put --acl-public ${PACKDIR}/${ARNAME} s3://rail-emr/bin/"
