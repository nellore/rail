#!/bin/sh

s3cmd put --acl-public bootstrap/* s3://tornado-emr/bootstrap/
