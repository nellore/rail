#!/bin/bash

# install-bowtie.sh
#
# Install SRA toolkit, including fastq-dump.  Version 2.1.7 as of now.

set -e

sudo apt-get --yes install sra-toolkit || { echo 'apt-get failed' ; exit 1; }
