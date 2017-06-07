#!/bin/sh

BUCKET=rail-emr

for i in *.zip ; do
    aws s3 cp "${i}" "s3://${BUCKET}/dependencies/" --acl public-read
done
