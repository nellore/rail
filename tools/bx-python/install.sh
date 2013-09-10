#!/bin/sh

hg clone http://bitbucket.org/james_taylor/bx-python/ bx-python-central
cd bx-python-central && sudo python setup.py install
