#!/usr/bin/env bash
# Make sure a minimal Elastic MapReduce cluster is running before using this script.
# $1 should be the Master public DNS of the cluster; e.g., ec2-54-92-221-31.compute-1.amazonaws.com
# $2 should be the full path to your pem file
# The reason for SSHing rather than doing this locally 
# is that we want to be sure we're using Amazon Hadoop's classes,
# including its particular hadoop-lzo.

CWD=$(pwd)
cd $(dirname "${BASH_SOURCE[0]}")
NEW_UUID=$(LC_CTYPE=C tr -dc A-Za-z < /dev/urandom | head -c 32 | xargs)
for JAR in relevant-elephant multiple-files mod-partitioner
do
	rm -f ${JAR}.tar.gz
	tar cvzf ${JAR}.tar.gz ${JAR}
	scp -i $2 ${JAR}.tar.gz hadoop@$1:~/${JAR}.tar.gz
done
ssh -t -t -i $2 hadoop@${1} <<ENDSSH
rm -f compile_${NEW_UUID}.log
for JAR in relevant-elephant multiple-files mod-partitioner
do
	tar xvzf \${JAR}.tar.gz
    rm -rf \${JAR}_out
    mkdir -p \${JAR}_out
    (javac -classpath \$(find ~/share/ *.jar | tr '\n' ':') -d \${JAR}_out \${JAR}/*.java 2>&1) >>compile_${NEW_UUID}.log
done
logout
ENDSSH
scp -i $2 hadoop@$1:~/compile_${NEW_UUID}.log ${CWD}/
for JAR in relevant-elephant multiple-files mod-partitioner
do
	rm -f ${JAR}.tar.gz
done
cd ${CWD}
echo 'Open' compile_${NEW_UUID}.log 'for compiler error logs.'