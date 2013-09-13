#!/usr/bin/perl

"""
reduce.py

Simple wrapper that mimics some of Hadoop's behavior during the Map step of a
MapReduce computation.
"""

import os
import sys
import argparse
import shutil
import time
import gzip
import bz2
import string
import multiprocessing
import subprocess
import heapq
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
    '--num-tasks', metavar='INT', type=int, required=True, help='Divide input into this many tasks.')
parser.add_argument(\
    '--bin-fields', metavar='INT', type=int, required=True, help='# fields to bin by.')
parser.add_argument(\
    '--sort-fields', metavar='INT', type=int, required=True, help='# fields to sort by.')
parser.add_argument(\
    '--external-sort', action='store_const', const=True, default=False, help='Use external program to sort tuples.')
parser.add_argument(\
    '--sort-size', metavar='INT', type=int, required=True, help='Memory cap on sorts.')
parser.add_argument(\
    '--messages', metavar='PATH', type=str, help='File to store stderr messages to')
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
    '--multiple-outputs', action='store_const', const=True, default=False, help='Outputs go in subdirectories labeled with first output token')
parser.add_argument(\
    '--keep-all', action='store_const', const=True, default=False, help='Keep all intermediate results')
parser.add_argument(\
    '--profile', action='store_const', const=True, default=False, help='Profile the code')
parser.add_argument(\
    '--verbose', action='store_const', const=True, default=False, help='Prints out extra debugging statements')

# Collect the reducer arguments
argv = sys.argv
reduceArgv = []
in_args = False
for i in xrange(1, len(sys.argv)):
    if in_args:
        reduceArgv.append(sys.argv[i])
    if sys.argv[i] == '--':
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
    intermediate = output + ".r.int"

message('==========================')
message('Step "%s" REDUCER' % args.name)
message('==========================')
message('Time: %s' % time.strftime('%H:%M:%S %d-%b-%Y'))
message('Inputs:')
for inp in args.input: message('  "%s"' % inp)
message('Output: "%s"' % output)
message('Intermediate: "%s"' % intermediate)
num_processes = args.num_processes or multiprocessing.cpu_count()
message('# parallel processes: %d' % num_processes)
message('Retries=%d, delay=%d seconds' % (args.num_retries, args.delay))
message('Bin fields=%d, sort fields=%d' % (args.bin_fields, args.sort_fields))
message('Max sort memory footprint=%d' % (args.sort_size))

options = []
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
    os.mkdir(d)
    if not os.path.exists(d) and os.path.isdir(d):
        mydie('Could not create new directory "%s"' % d, lev + 5)

# Where stdout output is dumped
checkDir(output, 100)

# Where stdout output is dumped
checkDir(intermediate, 200)

# Where stderr output is dumped
errDir = os.path.join(intermediate, 'reduce.err')
checkDir(errDir, 300)

# Where mapper programs are actually run
workingDir = os.path.join(intermediate, 'reduce.wds')
checkDir(workingDir, 400)

# Where stderr output from sort processes go
sortErrDir = os.path.join(intermediate, 'sort.err')
checkDir(sortErrDir, 500)

# Where unsorted reduce-task inputs are stored temporarily
taskDir = os.path.join(intermediate, 'reduce.tasks')
checkDir(taskDir, 600)

# Where sorted reduce task inputs are stored temporarily
staskDir = os.path.join(intermediate, 'reduce.stasks')
checkDir(staskDir, 700)

def openex(fn):
    if fn.endswith('.gz'): return gzip.GzipFile(fn, 'r')
    elif fn.endswith('.bz2'): return bz2.BZ2File(fn, 'r')
    else: return open(fn, 'r')

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

########################################
# Stage 1. Partition bins into tasks
########################################

# Go through all input files in parallel, calculating bin sizes for all bins.
binCountQueue = multiprocessing.Queue()
def binCountWorker(fn):
    cnt = defaultdict(int)
    with openex(fn) as fh:
        for ln in fh:
            ln = ln.rstrip()
            if len(ln) == 0: continue
            toks = string.split(ln, '\t')
            cnt['\t'.join(toks[:args.bin_fields])] += 1
    binCountQueue.put(cnt)

binCount = defaultdict(int)
countPool = multiprocessing.Pool(num_processes)
countPool.map(binCountWorker, inps)
checkFailQueue()

while not binCountQueue.empty():
    cnt = binCountQueue.get()
    for k, v in cnt.iteritems():
        binCount[k] += v

# Allocate bins to tasks, always adding to task with least tuples so far
taskNames = [ "task-%05d" % i for i in xrange(args.num_tasks) ]
taskq = [ (0, taskNames[i], []) for i in xrange(args.num_tasks) ]
keyToTask = {}
for k, v in binCount.iteritems():
    sz, nm, klist = heapq.heappop(taskq)
    sz += v
    klist.append(k)
    assert k not in keyToTask
    heapq.heappush(taskq, (sz, nm, klist))
    keyToTask[k] = nm

# Write out all the tasks to files within 'taskDir'
ofhs = {}
for inp in inps:
    with openex(inp) as fh:
        for ln in fh:
            ln = ln.rstrip()
            if len(ln) == 0: continue
            toks = string.split(ln, '\t')
            k = '\t'.join(toks[:args.bin_fields])
            task = keyToTask[k]
            if task not in ofhs:
                ofhs[task] = open(os.path.join(taskDir, task), 'w')
            ofhs[task].write(ln)
            ofhs[task].write('\n')

for fh in ofhs.itervalues():
    fh.close()

########################################
# Stage 2. Sort and reduce each task
########################################

def doSort(task, external=True, keep=False):
    assert external # only know how to use external sort for now
    inputFn, sortedFn = os.path.join(taskDir, task), os.path.join(staskDir, task)
    sortErrFn = os.path.join(sortErrDir, task)
    cmd = 'sort -S %d -k1,%d %s >%s 2>%s' % (args.sort_size, args.sort_fields, inputFn, sortedFn, sortErrFn)
    el = os.system(cmd)
    if el != 0:
        msg = 'Sort command "%s" for sort task "%s" failed with exitlevel: %d' % (cmd, task, el)
        failQ.put((msg, inputFn, sortErrFn, cmd))
    elif not keep:
        os.remove(inputFn)
        os.remove(sortErrFn)

sortPool = multiprocessing.Pool(num_processes)
sortPool.map(doSort, taskNames)
checkFailQueue()

reduceCmd = ' '.join(reduceArgv)

def doReduce(task, keep=False):
    sortedFn = os.path.join(staskDir, task)
    sortedFn = os.path.abspath(sortedFn)
    if not os.path.exists(sortedFn):
        raise RuntimeError('No such sorted task: "%s"' % sortedFn)
    outFn, errFn = os.path.join(output, task), os.path.join(errDir, task)
    outFn, errFn = os.path.abspath(outFn), os.path.abspath(errFn)
    cmd = 'cat %s | %s >%s 2>%s' % (sortedFn, reduceCmd, outFn, errFn)
    wd = os.path.join(workingDir, task)
    checkDir(wd, 800)
    pipe = subprocess.Popen(cmd, bufsize=-1, shell=True, cwd=wd)
    el = pipe.wait()
    if el != 0:
        msg = 'Reduce command "%s" for sort task "%s" failed with exitlevel: %d' % (cmd, task, el)
        failQ.put((msg, sortedFn, errFn, cmd))
    elif not keep:
        os.remove(sortedFn)
        os.remove(errFn)
        shutil.rmtree(wd)

reducePool = multiprocessing.Pool(num_processes)
reducePool.map(doReduce, taskNames)
checkFailQueue()

if not args.keep_all:
    message('Removing intermediate directory "%s"\n' % intermediate)
    shutil.rmtree(intermediate)

message('SUCCESS\n')
