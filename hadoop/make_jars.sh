#!/usr/bin/env bash
# Make sure a minimal Elastic MapReduce cluster is running before using this script.
# $1 should be the Master public DNS of the cluster; e.g., ec2-54-92-221-31.compute-1.amazonaws.com
# $2 should be the full path to your pem file
# The reason for SSHing rather than doing this locally 
# is that we want to be sure we're using Amazon Hadoop 2.4.0's classes,
# including its particular hadoop-lzo

cd $(dirname "${BASH_SOURCE[0]}")
rm -f relevantelephant.tar.gz relevantelephant.jar hadoop-streaming-mod.jar hadoop-streaming-mod.tar.gz
tar cvzf relevant-elephant.tar.gz relevant-elephant
tar cvzf hadoop-streaming-mod.tar.gz hadoop-streaming-mod
scp -i $2 relevant-elephant.tar.gz hadoop@$1:~/relevant-elephant.tar.gz
scp -i $2 hadoop-streaming-mod.tar.gz hadoop@$1:~/hadoop-streaming-mod.tar.gz
ssh -t -t -i $2 hadoop@${1} <<ENDSSH
tar xvzf relevant-elephant.tar.gz
rm -rf relevant-elephant_out
mkdir -p relevant-elephant_out
javac -classpath \$(find ~/share/ *.jar | tr '\n' ':') -d relevant-elephant_out relevant-elephant/*.java
jar -cvf relevant-elephant.jar -C relevant-elephant_out .
tar xvzf hadoop-streaming-mod.tar.gz
rm -rf hadoop-streaming-mod_out
mkdir -p hadoop-streaming-mod_out
javac -classpath \$(find ~/share/ *.jar | tr '\n' ':') -d hadoop-streaming-mod_out hadoop-streaming-mod/*.java
cp ~/share/hadoop/tools/lib/hadoop-streaming*.jar hadoop-streaming-mod_out/
cd hadoop-streaming-mod_out; mv hadoop-streaming*.jar hadoop-streaming-mod.jar; jar uf hadoop-streaming*.jar org/apache/hadoop/streaming/* edu/jhu/cs/MultipleOutputs*
mv hadoop-streaming-mod.jar ../
logout
ENDSSH
scp -i $2 hadoop@$1:~/hadoop-streaming-mod.jar ./
scp -i $2 hadoop@$1:~/relevant-elephant.jar ./
rm -rf ../lib
mkdir -p ../lib
cp hadoop-streaming-mod.jar ../lib/
cp relevant-elephant.jar ../lib/