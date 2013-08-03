#!/bin/bash

# install-python.sh
#
# Install Python and the key modules we use.

sudo apt-get --yes install python python-setuptools python-numpy python-scipy python-rpy2 python-argparse || { echo 'apt-get failed' ; exit 1; }
