#!/bin/bash

# install-kenttools.sh
#
# Install some of the more useful tools from Jim Kent for manipulating bed,
# wig, bigBed and bigWig files.

# Change to destination directory
mkdir -p $1
cd $1

if [ `uname -m` != "x86_64" ] ; then
	echo "Not a 64-bit platform!"
	exit 1
fi

# wget the needed files
for i in bedToBigBed bigBedToBed wigToBigWig bigWigToWig ; do
	wget -S -T 10 -t 5 http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/$i || { echo 'wget failed' ; exit 1; }
done
