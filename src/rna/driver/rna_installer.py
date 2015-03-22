#!/usr/bin/env python
"""
rna_installer.py
Part of Rail-RNA

Contains a class for installing Rail-RNA.
"""
import sys
from rna_config import *
base_path = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
import site
site.addsitedir(base_path)
import dependency_urls

class RailRnaInstaller(object):
    """ Installs Rail-RNA and its assorted dependencies. """
    def __init__(self):
        if sys.platform in ['linux', 'linux2']:
            bowtie1_url = http://downloads.sourceforge.net/project/bowtie-bio/bowtie/1.1.1/bowtie-1.1.1-linux-x86_64.zip
        elif sys.platform == 'darwin':

        else:
            print_to_screen(
                    'Rail-RNA cannot be installed because it is not supported '
                    'by your OS. Currently supported OSes are Mac OS X and '
                    'Linux.'
                )