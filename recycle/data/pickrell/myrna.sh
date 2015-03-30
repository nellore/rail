# ReCount-style Myrna command for analyzing Pickrell et al data, but with a
# bunch of spot instances

[ -z "$S3CFG" ] && echo "Set S3CFG" && exit 1
[ -z "$BUCKET" ] && echo "Set BUCKET" && exit 1
[ -z "$MYRNA_HOME" ] && echo "Set MYRNA_HOME" && exit 1

PROJ=pickrell
PROJ_FULL=Pickrell

s3cmd -c $S3CFG del --recursive s3://$BUCKET/myrna/$PROJ/output

$MYRNA_HOME/myrna_emr \
	--input=s3n://$BUCKET/myrna/$PROJ/preprocessed \
	--output=s3n://$BUCKET/myrna/$PROJ/output \
	--reference=s3n://myrna-refs/human_ensembl_61.jar \
	--instances=1,3,20 \
	--name "$PROJ_FULL" \
	--bowtie-args="-v 2 -m 1" \
	--bid-price=0.071 \
	--truncate=35 \
	--gene-footprint=intersect \
	--from-middle \
	$*
