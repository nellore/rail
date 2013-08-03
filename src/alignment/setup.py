#! /usr/bin/env python

# System imports
from distutils.core import *
from distutils      import sysconfig

# Third-party modules - we depend on numpy for everything
import numpy

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()
    
# inplace extension module
_nw = Extension("_nw",
                sources=["nw.i","nw.c"],
                include_dirs = [numpy_include],
                )

# NumyTypemapTests setup
setup(  name        = "needleman wunsch function",
        description = "performs an alignment between two sequences.",
        
        author      = "Jamie Morton",
        version     = "1.0",
        ext_modules = [_nw]
        )
