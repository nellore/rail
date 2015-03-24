# ReCount-style Myrna command for analyzing Geuvadis data with a
# bunch of spot instances

[ -z "$S3CFG" ] && echo "Set S3CFG" && exit 1
[ -z "$BUCKET" ] && echo "Set BUCKET" && exit 1
[ -z "$MYRNA_HOME" ] && echo "Set MYRNA_HOME" && exit 1

SMALL=
PROJ=geuvadis$SMALL
PROJ_FULL=Geuvadis$SMALL

s3cmd -c $S3CFG del --recursive s3://$BUCKET/myrna/$PROJ/output

# Cluster configuration I tried first was --instances=1,3,72
# Total of 75 worker/task nodes
# - 6 on-demand nodes, where price = $0.58 with $0.12 surcharge
# - 70 task nodes, where price = $0.071 with $0.12 surcharge
# - Total is ~$17.5 per hour for this cluster

$MYRNA_HOME/myrna_emr \
	--input=s3n://$BUCKET/myrna/$PROJ/preprocessed \
	--output=s3n://$BUCKET/myrna/$PROJ/output \
	--reference=s3n://myrna-refs/human_ensembl_61.jar \
	--instances=1,5,70 \
	--name "$PROJ_FULL" \
	--bowtie-args="-v 2 -m 1" \
	--bid-price=0.071 \
	--truncate=35 \
	--partition-len=3000 \
	--gene-footprint=intersect \
	--from-middle \
	--ditch-alignments \
	$*
