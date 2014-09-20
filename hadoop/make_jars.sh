#!/usr/bin/env bash
# Make sure a minimal Elastic MapReduce cluster is running before using this script.
# $1 should be the Master public DNS of the cluster; e.g., ec2-54-92-221-31.compute-1.amazonaws.com
# $2 should be the full path to your pem file
# The reason for SSHing rather than doing this locally 
# is that we want to be sure we're using Amazon Hadoop 2.4.0's classes,
# including its particular hadoop-lzo

cd $(dirname "${BASH_SOURCE[0]}")
rm -f relevantelephant.tar.gz relevantelephant.jar multiple-files.jar multiple-files.tar.gz
tar cvzf relevant-elephant.tar.gz relevant-elephant
tar cvzf multiple-files.tar.gz multiple-files
scp -i $2 relevant-elephant.tar.gz hadoop@$1:~/relevant-elephant.tar.gz
scp -i $2 multiple-files.tar.gz hadoop@$1:~/multiple-files.tar.gz
ssh -t -t -i $2 hadoop@${1} <<ENDSSH
tar xvzf relevant-elephant.tar.gz
rm -rf relevant-elephant_out
mkdir -p relevant-elephant_out
javac -classpath \$(find ~/share/ *.jar | tr '\n' ':') -d relevant-elephant_out relevant-elephant/*.java
jar -cvf relevant-elephant.jar -C relevant-elephant_out .
tar xvzf multiple-files.tar.gz
rm -rf multiple-files_out
mkdir -p multiple-files_out
javac -classpath \$(find ~/share/ *.jar | tr '\n' ':') -d multiple-files_out multiple-files/*.java
jar -cvf multiple-files.jar -C multiple-files_out .
logout
ENDSSH
scp -i $2 hadoop@$1:~/multiple-files.jar ./
scp -i $2 hadoop@$1:~/relevant-elephant.jar ./
rm -rf ../lib
mkdir -p ../lib
cp multiple-files.jar ../lib/
cp relevant-elephant.jar ../lib/