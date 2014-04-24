#!/usr/bin/env python
"""
bowtie.py

Alignment module 
"""

import os
import sys
import subprocess
import threading
import site

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, 'util'))

import path

def addArgs(parser):
    parser.add_argument(\
        '--bowtie-exe', metavar='EXE', type=str, required=False,
        default='bowtie',
        help='Path to executable for Bowtie')
    parser.add_argument(\
        '--bowtie-build-exe', metavar='EXE', type=str, required=False,
        default='bowtie-build',
        help='Path to executable for Bowtie-build')
    parser.add_argument(\
        '--bowtie-idx', metavar='INDEX', type=str, required=False,
        default='',
        help='Path to Bowtie index. Specify its basename.')
    parser.add_argument(\
        '--bowtie2-exe', metavar='EXE', type=str, required=False,
        default='bowtie2',
        help='Path to executable for Bowtie2')
    parser.add_argument(\
        '--bowtie2-build-exe', metavar='EXE', type=str, required=False,
        default='bowtie2-build',
        help='Path to executable for Bowtie2-build')
    parser.add_argument(\
        '--bowtie2-idx', metavar='INDEX', type=str, required=False,
        default='',
        help='Path to Bowtie2 index. Specify its basename.')

def out(pi):
    for line in iter(pi.readline, ''):
        print line

def cmd(bowtieExe="bowtie", bowtieIdx="genome", readFn=None, bowtieArgs=None,
        sam=False, outputFilename=''):
    # Check that Bowtie exists and is executable
    if not path.is_exe(bowtieExe) and path.which(bowtieExe) is None:
        print >>sys.stderr, "Bowtie executable '%s' cannot be run" % bowtieExe
        sys.exit(1)
    # Check that index is there
    for i in [".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt",
                ".rev.1.ebwt", ".rev.2.ebwt"]:
        if not os.path.isfile(bowtieIdx + i):
            print >>sys.stderr, "Could not find index file " + (bowtieIdx + i)
            sys.exit(1)
    out_arg = ""
    if sam:
        out_arg = " -S "
    argstr = ''
    mycmd = "%s %s --mm %s %s --12 " % (bowtieExe, bowtieArgs,
                                            out_arg, bowtieIdx)
    if readFn is None:
        mycmd += '-'
    else:
        mycmd += readFn + ' ' + outputFilename
    return mycmd

def proc(bowtieExe="bowtie", bowtieIdx="genome", readFn=None, bowtieArgs=None,
    sam=False, stdoutPipe=False, outHandler=None, errHandler=None,
    stdinPipe=True, outputFilename=''):
    stdout_pipe = None if (outHandler is None
                                and not stdoutPipe) else subprocess.PIPE
    stderr_pipe = None if errHandler is None else subprocess.PIPE
    stdin_pipe = subprocess.PIPE if stdinPipe else None
    mycmd = cmd(bowtieExe=bowtieExe, bowtieIdx=bowtieIdx, readFn=readFn,
                bowtieArgs=bowtieArgs, sam=sam, outputFilename=outputFilename)
    print >> sys.stderr, "Starting command: '%s'" % mycmd
    proc = subprocess.Popen(\
        mycmd, shell=True, stdin=stdin_pipe, stdout=stdout_pipe,
        stderr=stderr_pipe, bufsize=-1)
    threads = []
    if outHandler is not None:
        print >> sys.stderr, "  Starting stdout handler"
        t = threading.Thread(target=outHandler, args=(proc.stdout,))
        t.daemon = True # thread dies with the program
        t.start()
        threads.append(t)
    if errHandler is not None:
        print >> sys.stderr, "  Starting stderr handler"
        t = threading.Thread(target=errHandler, args=(proc.stderr,))
        t.daemon = True # thread dies with the program
        t.start()
        threads.append(t)
    return proc, mycmd, threads
