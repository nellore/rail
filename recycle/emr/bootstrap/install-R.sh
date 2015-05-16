#!/bin/bash

# install-R.sh
#
# Install r-base.  Seems to be version 2.11.1.

set -e

sudo apt-get --yes install r-base || { echo 'apt-get failed' ; exit 1; }
