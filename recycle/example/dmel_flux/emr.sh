#!/bin/sh

BUCKET=langmead

if [ "$1" = "upload" ] ; then
	s3cmd put *.fastq s3://$BUCKET/dmel_flux/fastq/
	
	cat >dmel_flux.s3.manifest <<EOF
s3://$BUCKET/dmel_flux/fastq/dmel_flux-0.fastq	0	dmel_flux-0-0
s3://$BUCKET/dmel_flux/fastq/dmel_flux-1.fastq	0	dmel_flux-0-1
s3://$BUCKET/dmel_flux/fastq/dmel_flux-2.fastq	0	dmel_flux-0-2
s3://$BUCKET/dmel_flux/fastq/dmel_flux-3.fastq	0	dmel_flux-0-3
s3://$BUCKET/dmel_flux/fastq/dmel_flux-4.fastq	0	dmel_flux-0-4
s3://$BUCKET/dmel_flux/fastq/dmel_flux-5.fastq	0	dmel_flux-0-5
EOF
	s3cmd put dmel_flux.s3.manifest s3://$BUCKET/dmel_flux/manifest/
fi

s3cmd del --recursive s3://$BUCKET/dmel_flux/intermediate
s3cmd del --recursive s3://$BUCKET/dmel_flux/output

python $RAIL_HOME/src/driver/rail-rna.py \
	--emr \
	--no-differential \
	--manifest s3://$BUCKET/dmel_flux/manifest/dmel_flux.s3.manifest \
	--input s3://$BUCKET/dmel_flux/fastq \
	--intermediate s3://$BUCKET/dmel_flux/intermediate \
	--output s3://$BUCKET/dmel_flux/output \
	--reference s3://tornado-emr/refs/dm3_UCSC.tar.gz \
	--instance-type c1.xlarge \
	--instance-counts 1,0,0 \
	$*
