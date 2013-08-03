#!/bin/bash

# install-pypy.sh
#
# Install PyPy.

# Got this list from http://doc.pypy.org/en/latest/getting-started-python.html
sudo apt-get --yes install \
 gcc make python-dev libffi-dev libsqlite3-dev pkg-config \
 libz-dev libbz2-dev libncurses-dev libexpat1-dev \
 libssl-dev libgc-dev python-sphinx python-greenlet || { echo "apt-get failed" ; exit 1 }

# TODO: actually get PyPy.  Note that the pypy-c binary is not enough, you
# also need to have all the supporting directories!
