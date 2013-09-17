#!/usr/bin/perl

"""
map.py

Simple wrapper that mimics some of Hadoop's behavior during the Map step of a
MapReduce computation.
"""

import os
import sys
import argparse
import shutil
import time
import subprocess
import gzip
import bz2
import string
import multiprocessing
from collections import defaultdict

parser = argparse.ArgumentParser(description=\
    'Simple wrapper that mimics some of Hadoop\'s behavior during the Map step of a MapReduce computation.')

parser.add_argument(\
    '--name', metavar='PATH', type=str, required=True, help='Name of this step')
parser.add_argument(\
    '--input', metavar='PATH', type=str, required=True, help='Files with input data', nargs='+')
parser.add_argument(\
    '--output', metavar='PATH', type=str, required=True, help='Files with output data')
parser.add_argument(\
    '--messages', metavar='PATH', type=str, help='File to store stderr messages to')
parser.add_argument(\
    '--counters', metavar='PATH', type=str, help='File to store counter data to')
parser.add_argument(\
    '--intermediate', metavar='PATH', type=str, help='Directory to store intermediate data in')
parser.add_argument(\
    '--num-processes', metavar='INT', type=int, help='Max # of simultaneous processes to run.  Default = # processors you have.')
parser.add_argument(\
    '--num-retries', metavar='INT', type=int, default=3, help='# times to retry the mapper if something goes wrong')
parser.add_argument(\
    '--delay', metavar='INT', type=int, default=5, help='# seconds to wait before a retry')
parser.add_argument(\
    '--force', action='store_const', const=True, default=False, help='Profile the code')
parser.add_argument(\
    '--line-by-line', action='store_const', const=True, default=False, help='Process each line as a new')
parser.add_argument(\
    '--multiple-outputs', action='store_const', const=True, default=False, help='Outputs go in subdirectories labeled with first output token')
parser.add_argument(\
    '--keep-all', action='store_const', const=True, default=False, help='Keep all intermediate results')
parser.add_argument(\
    '--profile', action='store_const', const=True, default=False, help='Profile the code')
parser.add_argument(\
    '--verbose', action='store_const', const=True, default=False, help='Prints out extra debugging statements')

# Collect the bowtie arguments first
argv = sys.argv[:]
mapargv = []
in_args = False
for i in xrange(1, len(sys.argv)):
    if in_args:
        mapargv.append(sys.argv[i])
    elif sys.argv[i] == '--':
        argv = sys.argv[:i]
        in_args = True

args = parser.parse_args(argv[1:])

msgfhs = [ sys.stderr ]

def message(s, add_nl=True):
    for fh in msgfhs:
        fh.write(s)
        if add_nl: fh.write('\n')

def mydie(msg, lev):
    message("Fatal error %d:\n%s" % (lev, msg))
    sys.exit(lev)

if args.messages: msgfhs.append(open(args.messages, 'w'))

inps = args.input
for inp in inps:
    if not os.path.exists(inp):
        raise RuntimeError("--input doesn't exist: \"%s\"" % inp)
output = args.output
intermediate = args.intermediate
if intermediate is None:
    intermediate = output + ".m.int"

message('==========================')
message('Step "%s" MAPPER' % args.name)
message('==========================')
message('Time: %s' % time.strftime('%H:%M:%S %d-%b-%Y'))
message('Inputs:')
for inp in args.input: message('  "%s"' % inp)
message('Output: "%s"' % output)
message('Intermediate: "%s"' % intermediate)
num_processes = args.num_processes or multiprocessing.cpu_count()
assert num_processes > 0
message('# parallel processes: %d' % num_processes)
message('Retries=%d, delay=%d seconds' % (args.num_retries, args.delay))

options = []
if args.line_by_line: options.append("--line-by-line")
if args.keep_all: options.append("--keep-all")
if args.force: options.append("--force")
if args.multiple_outputs: options.append("--multiple-outputs")
message('Options: [' + ' '.join(options) + ']')

def checkDir(d, lev):
    """ Check whether a directory exists.  If so and --force is specified,
        remove it.  Create it if it doesn't exist. """
    if os.path.exists(d):
        if args.force:
            message('Removing "%s" due to --force' % d)
            shutil.rmtree(d)
        else:
            mydie('Output directory "%s" already exists' % d, lev)
    os.makedirs(d)
    if not os.path.exists(d) and os.path.isdir(d):
        mydie('Could not create new directory "%s"' % d, lev + 5)

# Where stdout output is dumped
checkDir(output, 100)

# Where stdout output is dumped
checkDir(intermediate, 200)

# Where stderr output is dumped
errDir = os.path.join(intermediate, 'map.err')
checkDir(errDir, 300)

# Where mapper programs are actually run
workingDir = os.path.join(intermediate, 'map.wds')
checkDir(workingDir, 400)

countfh = None
if args.counters: countfh = open(args.counters, 'w')

def openex(fn):
    if fn.endswith('.gz'): return gzip.GzipFile(fn, 'r')
    elif fn.endswith('.bz2'): return bz2.BZ2File(fn, 'r')
    else: return open(fn, 'r')

def fileizeInput(inps):
    """ If any inputs are directories, replace the directory with all the
        files within. """
    newinps = []
    for inp in inps:
        if os.path.isdir(inp):
            for fn in os.listdir(inp):
                fn = os.path.join(inp, fn)
                if os.path.isfile(fn):
                    newinps.append(fn)
        else:
            newinps.append(inp)
    return newinps

assert len(inps) > 0
inps = fileizeInput(inps)
assert len(inps) > 0

# If input is line-by-line, read lines into lineinp list
lineinp = []
for inp in map(os.path.abspath, inps):
    if not os.path.exists(inp):
        mydie('No such input file "%s"' % inp, 500)
    if not os.path.isfile(inp):
        mydie('Input "%s" is not a file' % inp, 600)
    if args.line_by_line:
        with openex(inp) as ifh:
            for ln in ifh:
                ln = ln.rstrip()
                if len(ln) == 0 or ln[0] == '#':
                    continue
                lineinp.append(ln)

def mkdirQuiet(d, lev):
    try: os.makedirs(d)
    except OSError: pass
    if not os.path.exists(d) and os.path.isdir(d):
        mydie('Could not create directory "%s"' % d, lev)

failQ = multiprocessing.Queue()

def checkFailQueue():
    if not failQ.empty():
        while not failQ.empty():
            err, inp, errfn, cmd = failQ.get()
            message('******')
            message('* ' + err)
            message('* Command was:')
            message('*   ' + cmd)
            message('* Input file/string was:')
            message('*   ' + inp)
            message('* Error message is in file:')
            message('*   ' + errfn)
            # TODO: echo error message from error file
            message('******')
        message('FAILED')
        sys.exit(1)

cmd = ' '.join(mapargv)
taskn = len(inps)

def doMapper(tup):
    """ Run a single mapper task """
    inp, taski = tup
    if not failQ.empty():
        return None
    message('Pid %d processing input "%s" [%d of %d]' % (os.getpid(), inp, taski, taskn))
    ofn = "map-%05d" % taski
    mkdirQuiet(output, 700)
    mkdirQuiet(errDir, 800)
    outFullFn = os.path.join(output, ofn)
    errFullFn = os.path.join(errDir, ofn)
    mycmd = cmd + " >%s 2>%s" % (outFullFn, errFullFn)
    wd = os.path.join(workingDir, str(taski))
    mkdirQuiet(wd, 900) # make the working directory
    ret = 0
    fullcmd = mycmd
    if not args.line_by_line:
        # TODO: perhaps handle gzip or bzip2 internally in Python instead of
        # setting up a pipeline of external tools like this
        if inp.endswith('.gz'):
            fullcmd = ('gzip -dc %s | ' % inp) + mycmd
        elif inp.endswith('.bz2'):
            fullcmd = ('bzip2 -dc %s | ' % inp) + mycmd
        else:
            fullcmd = ('cat %s | ' % inp) + mycmd
    for _ in xrange(args.num_retries+1):
        if args.line_by_line:
            pipe = subprocess.Popen(mycmd, bufsize=-1, stdin=subprocess.PIPE, shell=True, cwd=wd)
            pipe.stdin.write(inp + '\n')
            pipe.stdin.close()
            ret = pipe.wait()
            if ret == 0:
                return outFullFn
            message('Non-zero return (%d) after closing pipe "%s"' % (ret, mycmd))
        else:
            ret = os.system(fullcmd)
            if ret == 0:
                return outFullFn
            message('Non-zero return (%d) after executing "%s"' % (ret, fullcmd))
        message('Retrying in %d seconds...' % args.delay)
        time.sleep(args.delay)
    # Finished
    failQ.put([\
        "Mapper %d of %d (pid %d) failed the maximum # of times %d" % (taski, taskn, os.getpid(), args.num_retries+1),
        inp, errFullFn, fullcmd])
    return outFullFn

num_processes = min(num_processes, len(inps))
message('Starting %d processes with command: "%s"' % (num_processes, cmd))

tasks = []
taski, taskn = 1, len(inps)
for inp in map(os.path.abspath, inps):
    tasks.append((inp, taski))
    taski += 1

pool = multiprocessing.Pool(num_processes)
outfns = pool.map(doMapper, tasks)
checkFailQueue()

if args.multiple_outputs:
    
    def doSplit(tup, keep=args.keep_all):
        task, outFn = tup
        splitDir = os.path.join(output, '_'.join([task, 'split']))
        checkDir(splitDir, 900)
        kToFh, kToFn = {}, {}
        with openex(outFn) as fh:
            for ln in fh:
                ln = ln.rstrip()
                if len(ln) == 0: continue
                toks = string.split(ln, '\t')
                if toks[0] not in kToFh:
                    fn = os.path.join(splitDir, '_'.join([toks[0], task]))
                    kToFn[toks[0]] = fn
                    kToFh[toks[0]] = open(fn, 'w')
                kToFh[toks[0]].write('\t'.join(toks[1:]))
                kToFh[toks[0]].write('\n')
        for ofh in kToFh.itervalues():
            ofh.close()
        if not keep:
            os.remove(outFn)
        return kToFn
    
    splitTups = []
    i = 0
    for outfn in outfns:
        splitTups.append(('task-%05d' % i, outfn))
        i += 1
    splitPool = multiprocessing.Pool(num_processes)
    kToFns = splitPool.map(doSplit, splitTups)
    checkFailQueue()
    kToFnList = defaultdict(list)
    for kToFn in kToFns:
        for k, v in kToFn.iteritems():
            kToFnList[k].append(v)
    
    # Now join them back up
    for k, vl in kToFnList.iteritems():
        kOutDir = os.path.join(output, k)
        checkDir(kOutDir, 900)
        for v in vl:
            fn = os.path.basename(v)
            shutil.copyfile(v, os.path.join(kOutDir, fn))
            if not args.keep_all:
                os.remove(v)
    
    # Remove all the split directories
    if not args.keep_all:
        removedAlready = set()
        for k, vl in kToFnList.iteritems():
            for v in vl:
                vdir = os.path.dirname(v)
                if vdir not in removedAlready:
                    shutil.rmtree(vdir)
                    removedAlready.add(vdir)

if not args.keep_all:
    message('Removing intermediate directory "%s"\n' % intermediate)
    shutil.rmtree(intermediate)

message('SUCCESS\n')
