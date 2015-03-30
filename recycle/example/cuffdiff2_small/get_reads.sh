#!/bin/sh

mkdir -p reads
cd reads
s3cmd $* get s3://tornado-emr/reads/cuffdiff2_small/* .
cd ..
