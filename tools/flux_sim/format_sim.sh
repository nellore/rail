#!/bin/bash
#appends the file name in front of every read name
#ex: sh format_sim.sh small_test

directory=$1
for FASTQ in $directory/*.fastq
do
    name=`basename $FASTQ .fastq`
    echo $name
    sed -i '1~4 s/\@/\@'$name'./' $FASTQ
done