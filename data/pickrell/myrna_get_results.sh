[ -z "$S3CFG" ] && echo "Set S3CFG" && exit 1
[ -z "$BUCKET" ] && echo "Set BUCKET" && exit 1

s3cmd -c $S3CFG get s3://$BUCKET/myrna/pickrell/output/myrna_results/results.tar.gz
