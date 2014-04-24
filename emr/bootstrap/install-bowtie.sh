#!/bin/bash

# install-bowtie.sh
#
# Install bowtie version 1.0.1 as of now.

set -e

sudo apt-get --yes install bowtie || { echo 'apt-get failed' ; exit 1; }
