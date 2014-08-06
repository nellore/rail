#!/bin/bash

# install-bowtie2.sh
#
# Install bowtie version 2.2.2 as of now.

set -e

wget http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.2.3/bowtie2-2.2.3-linux-x86_64.zip || { echo 'wget failed' ; exit 1; }
unzip bowtie2-2.2.3-linux-x86_64.zip || { echo 'unzip failed' ; exit 1; }
sudo ln -s /home/hadoop/bowtie2-2.2.3/bowtie2 /bin/bowtie2
sudo ln -s /home/hadoop/bowtie2-2.2.3/bowtie2-build /bin/bowtie2-build