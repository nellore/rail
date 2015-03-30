#!/bin/sh

[ -z "$TORNADO_HOME" ] && echo "Must set TORNADO_HOME" && exit 1

cat MANIFEST.myrna | \
python ../../src/rnawesome/preprocess.py \
	--gzip-output \
	--s3cfg $HOME/.ec2_jhu/.s3cfg \
	--acl-public \
	--push s3://tornado-emr/reads/bodymap2/ \
	$*
