#!/bin/sh

# build_pypy.sh
#
# Download pypy source and build a binary.  This is presumably something we do
# on an EMR machine so that there's a working version built on the appropriate
# OS.
#
# BTL: This took a very long time when I last tried it, and I'm not confident
# that what gets produced is portable.  If anyone discovers a good way to get
# pypy up and running on a new EMR instance quickly (such that we could do so
# in a bootstrap action) please let me know.

# Got this list from http://doc.pypy.org/en/latest/getting-started-python.html
sudo apt-get --yes install \
 gcc make python-dev libffi-dev libsqlite3-dev pkg-config \
 libz-dev libbz2-dev libncurses-dev libexpat1-dev \
 libssl-dev libgc-dev python-sphinx python-greenlet || { echo "apt-get failed" ; exit 1 }

VER=3.1
NM=pypy-${VER}-src
AR=${NM}.tar.bz2

wget https://bitbucket.org/pypy/pypy/downloads/${AR}

tar zxvj ${AR}
[ -d ${NM} ] || { echo "Directory ${NM} does not exist" ; exit 1 }
cd ${NM}/pypy/goal
python ../../rpython/bin/rpython --opt=jit targetpypystandalone.py || { echo "pypy build failed" ; exit 1}
