#!/usr/bin/env python

"""
align.py

Alignment module 
"""

import os
import sys
import subprocess
import threading

def addArgs(parser):
    parser.add_argument(\
        '--bowtieArgs', metavar='ARGS', type=str, nargs='+',
         required=False, help='Arguments to pass to Bowtie')
    parser.add_argument(\
        '--bowtieExe', metavar='EXE', type=str, required=False,
        help='Path to executable for Bowtie')
    parser.add_argument(\
        '--bowtieIdx', metavar='INDEX', type=str, required=False,
        help='Path to Bowtie index.  Specify its basename.')

def out(pi):
    for line in iter(pi.readline, ''):
        print line

def cmd(args, sam=False):
    # Check that Bowtie exists and is executable
    if not os.path.isfile(args.bowtieExe):
        print >>sys.stderr, "Bowtie executable '%s' does not exist" % args.bowtieExe
        sys.exit(1)
    if not os.access(args.bowtieExe, os.X_OK):
        print >>sys.stderr, "Bowtie executable '%s' exists but is not executable" % args.bowtieExe
        sys.exit(1)
    # Check that index is there
    for i in [".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt"]:
        if not os.path.isfile(args.bowtieIdx + i):
            print >>sys.stderr, "Could not find index file " + (args.bowtieIdx + i)
            sys.exit(1)
    out_arg = ""
    if sam:
        out_arg = " -S "
    return "%s %s --mm %s %s --12 -" % (args.bowtieExe, ' '.join(args.bowtieArgs), out_arg, args.bowtieIdx)

def proc(args, sam=False, outHandler=None, errHandler=None):
    stdout_pipe = None if outHandler is None else subprocess.PIPE
    stderr_pipe = None if errHandler is None else subprocess.PIPE
    proc = subprocess.Popen(\
        cmd(args, sam=sam),
        shell=True, stdin=subprocess.PIPE, stdout=stdout_pipe, stderr=stderr_pipe)
    if outHandler is not None:
        t = threading.Thread(target=outHandler, args=(proc.stdout,))
        t.daemon = True # thread dies with the program
        t.start()
    if errHandler is not None:
        t = threading.Thread(target=errHandler, args=(proc.stderr,))
        t.daemon = True # thread dies with the program
        t.start()
    return proc
