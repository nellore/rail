#!/bin/bash

# install-python.sh
#
# Install Python and the key modules we use.

sudo apt-get --yes install python-setuptools || { echo 'apt-get failed' ; exit 1; }
for mid in numpy scipy argparse rpy2 ; do
	sudo easy_install $mid || { echo "easy_install $mid failed" ; exit 1; }
done
