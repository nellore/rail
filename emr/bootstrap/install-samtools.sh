#!/bin/bash

# install-samtools.sh
#
# Install samtools version 0.1.8-1 as of now.

set -e

# samtools requires curses.h
sudo yum -y install ncurses-devel ncurses || { echo 'curses installation failed' ; exit 1; }
wget http://downloads.sourceforge.net/project/samtools/samtools/0.1.19/samtools-0.1.19.tar.bz2 || { echo 'wget failed' ; exit 1; }
tar xvjf samtools-0.1.19.tar.bz2 || { echo 'unzip failed' ; exit 1; }
cd samtools-0.1.19
sudo make || { echo 'make failed'; exit 1; }
cd ..
sudo ln -s /home/hadoop/samtools-0.1.19/samtools /bin/samtools