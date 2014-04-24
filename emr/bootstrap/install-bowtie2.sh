#!/bin/bash

# install-bowtie2.sh
#
# Install bowtie version 2.2.2 as of now.

set -e

sudo apt-get --yes install bowtie2 || { echo 'apt-get failed' ; exit 1; }
