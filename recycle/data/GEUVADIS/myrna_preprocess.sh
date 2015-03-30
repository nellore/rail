#!/bin/sh

[ -z "$S3CFG" ] && echo "Set S3CFG" && exit 1
[ -z "$BUCKET" ] && echo "Set BUCKET" && exit 1
[ -z "$MYRNA_HOME" ] && echo "Set MYRNA_HOME" && exit 1

PROJ=geuvadis
PROJ_FULL=Geuvadis
MANIFEST=manifest_mate1_randomized

[ ! -f "$MANIFEST" ] && echo "Could not find manifest file $MANIFEST" && exit 1

s3cmd del --recursive s3://$BUCKET/myrna/${PROJ}/preprocessed
s3cmd put $MANIFEST s3://$BUCKET/myrna/${PROJ}/manifest/

perl $MYRNA_HOME/myrna_emr \
	--input=s3n://$BUCKET/myrna/${PROJ}/manifest/$MANIFEST \
	--output=s3n://$BUCKET/myrna/${PROJ}/preprocessed \
	--just-preprocess \
	--instances=1,1,5 \
	--instance-type c1.xlarge \
	--name "${PROJ_FULL}Preprocess" \
	--bid-price=0.071 \
	$*
