#!/usr/bin/env python

"""
align.py

Alignment module 
"""

import os
import sys
import subprocess
import threading
import site

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "util"))

import path

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

def cmd(args, readFn=None, bowtieArgs=None, sam=False):
    # Check that Bowtie exists and is executable
    bowtieExe = args.bowtieExe or "bowtie"
    if not path.is_exe(bowtieExe) and path.which(bowtieExe) is None:
        print >>sys.stderr, "Bowtie executable '%s' cannot be run" % bowtieExe
        sys.exit(1)
    # Check that index is there
    for i in [".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt"]:
        if not os.path.isfile(args.bowtieIdx + i):
            print >>sys.stderr, "Could not find index file " + (args.bowtieIdx + i)
            sys.exit(1)
    out_arg = ""
    if sam:
        out_arg = " -S "
    argstr = ''
    if bowtieArgs is not None:
        argstr = ' '.join(bowtieArgs)
    if args.bowtieArgs is not None:
        argstr += ' '
        argstr += ' '.join(args.bowtieArgs)
    mycmd = "%s %s --mm %s %s --12 " % (bowtieExe, argstr, out_arg, args.bowtieIdx)
    mycmd += ("-" if readFn is None else readFn)
    return mycmd

def proc(args, readFn=None, bowtieArgs=None, sam=False, outHandler=None, errHandler=None, stdinPipe=True):
    stdout_pipe = None if outHandler is None else subprocess.PIPE
    stderr_pipe = None if errHandler is None else subprocess.PIPE
    stdin_pipe = subprocess.PIPE if stdinPipe else None
    mycmd = cmd(args, readFn=readFn, bowtieArgs=bowtieArgs, sam=sam)
    print >> sys.stderr, "Starting command: '%s'" % mycmd
    proc = subprocess.Popen(\
        mycmd, shell=True, stdin=stdin_pipe, stdout=stdout_pipe, stderr=stderr_pipe, bufsize=-1)
    if outHandler is not None:
        print >> sys.stderr, "  Starting stdout handler"
        t = threading.Thread(target=outHandler, args=(proc.stdout,))
        t.daemon = True # thread dies with the program
        t.start()
    if errHandler is not None:
        print >> sys.stderr, "  Starting stderr handler"
        t = threading.Thread(target=errHandler, args=(proc.stderr,))
        t.daemon = True # thread dies with the program
        t.start()
    return proc, mycmd
