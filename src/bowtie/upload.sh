#!/bin/sh

for i in *.zip ; do
    aws s3 cp "${i}" s3://rail-emr-requester-pays/dependencies/ --acl public-read
done
