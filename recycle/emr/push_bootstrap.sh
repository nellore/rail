#!/bin/sh

s3cmd put --acl-public bootstrap/* s3://rail-emr/bootstrap/
