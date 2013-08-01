#!/bin/sh

if [ -z "$TORNADO_HOME" ] ; then
	echo "Set TORNADO_HOME first"
	exit 1
fi

s3cmd del --recursive s3://langmead/tornado_cuffdiff2

python $TORNADO_HOME/src/driver/tornado.py \
	--emr \
	--input s3://tornado-emr/reads/cuffdiff2 \
	--output s3://langmead/tornado_cuffdiff2 \
	--reference s3://tornado-emr/refs/hg19_UCSC.tar.gz \
	--just-align \
	--instance-type c1.xlarge \
	--instance-counts 1,3,0 \
	$*

echo "s3cmd del --recursive s3://langmead/tornado_cuffdiff2"
