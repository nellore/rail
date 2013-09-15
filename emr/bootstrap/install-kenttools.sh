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
	# Note: the binaries in the linux.x86_64 directory stopped working as of
	# summer 2013.  I hope the linux.x86_64.v287 directory sticks around
	# forever!
	wget -S -T 10 -t 5 http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v287/$i || { echo 'wget failed' ; exit 1; }
	chmod a+x $i
	if [ ! -x $i ] ; then
		echo "$1/$i does not exist or isn't executable"
		exit 1
	fi
done
