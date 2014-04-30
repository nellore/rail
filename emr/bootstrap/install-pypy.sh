#!/bin/bash

# install-pypy.sh
#
# Install PyPy as of now.
# Reference: http://stackoverflow.com/questions/21803988/pypy-issue-with-shared-libraries-libffi-so-5

set -e
sudo add-apt-repository --yes ppa:pypy/ppa
sudo apt-get --yes install pypy pypy || { echo 'apt-get failed'; exit 1; }