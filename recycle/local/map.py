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

parser = argparse.ArgumentParser(
    description='Simple wrapper that mimics some of Hadoop\'s behavior during the Map step of a MapReduce computation.')

parser.add_argument('--name', metavar='PATH', type=str, required=True, help='Name of this step')
parser.add_argument('--input', metavar='PATH', type=str, required=True, help='Files with input data', nargs='+')
parser.add_argument('--output', metavar='PATH', type=str, required=True, help='Files with output data')
parser.add_argument('--messages', metavar='PATH', type=str, help='File to store stderr messages to')
parser.add_argument('--counters', metavar='PATH', type=str, help='File to store counter data to')
parser.add_argument('--intermediate', metavar='PATH', type=str, help='Directory to store intermediate data in')
parser.add_argument('--num-processes', metavar='INT', type=int,
                    help='Max # of simultaneous processes to run.  Default = # processors you have.')
parser.add_argument('--num-retries', metavar='INT', type=int, default=3,
                    help='# times to retry the mapper if something goes wrong')
parser.add_argument('--delay', metavar='INT', type=int, default=5, help='# seconds to wait before a retry')
parser.add_argument('--force', action='store_const', const=True, default=False, help='Profile the code')
parser.add_argument('--line-by-line', action='store_const', const=True, default=False,
                    help='Process each line as a new')
parser.add_argument('--multiple-outputs', action='store_const', const=True, default=False,
                    help='Outputs go in subdirectories labeled with first output token')
parser.add_argument('--keep-all', action='store_const', const=True, default=False, help='Keep all intermediate results')
parser.add_argument('--profile', action='store_const', const=True, default=False, help='Profile the code')
parser.add_argument('--verbose', action='store_const', const=True, default=False,
                    help='Prints out extra debugging statements')

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

msgfhs = [sys.stderr]


def message(s, add_nl=True):
    for fh in msgfhs:
        fh.write(s)
        if add_nl:
            fh.write('\n')


def mydie(msg, lev):
    message("Fatal error %d:\n%s" % (lev, msg))
    sys.exit(lev)


if args.messages:
    msgfhs.append(open(args.messages, 'w'))

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
for inp in args.input:
    message('  "%s"' % inp)
message('Output: "%s"' % output)
message('Intermediate: "%s"' % intermediate)
num_processes = args.num_processes or multiprocessing.cpu_count()
assert num_processes > 0
message('# parallel processes: %d' % num_processes)
message('Retries=%d, delay=%d seconds' % (args.num_retries, args.delay))

options = []
if args.line_by_line:
    options.append("--line-by-line")
if args.keep_all:
    options.append("--keep-all")
if args.force:
    options.append("--force")
if args.multiple_outputs:
    options.append("--multiple-outputs")
message('Options: [' + ' '.join(options) + ']')


def check_dir(d, lev):
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
check_dir(output, 100)

# Where stdout output is dumped
check_dir(intermediate, 200)

# Where stderr output is dumped
errDir = os.path.join(intermediate, 'map.err')
check_dir(errDir, 300)

# Where mapper programs are actually run
workingDir = os.path.join(intermediate, 'map.wds')
check_dir(workingDir, 400)

countfh = None
if args.counters:
    countfh = open(args.counters, 'w')


def openex(fn):
    if fn.endswith('.gz'):
        return gzip.GzipFile(fn, 'r')
    elif fn.endswith('.bz2'):
        return bz2.BZ2File(fn, 'r')
    else:
        return open(fn, 'r')


def fileize_input(inps):
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
inps = fileize_input(inps)
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


def mkdir_quiet(d, lev):
    try:
        os.makedirs(d)
    except OSError:
        pass
    if not os.path.exists(d) and os.path.isdir(d):
        mydie('Could not create directory "%s"' % d, lev)


failQ = multiprocessing.Queue()


def check_fail_queue():
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


def do_mapper(tup):
    """ Run a single mapper task """
    name, inp, taski = tup
    if not failQ.empty():
        return None
    message('Pid %d processing input "%s" [%d of %d]' % (os.getpid(), name, taski, taskn))
    ofn = "map-%05d" % taski
    mkdir_quiet(output, 700)
    mkdir_quiet(errDir, 800)
    out_full_fn = os.path.abspath(os.path.join(output, ofn))
    err_full_fn = os.path.abspath(os.path.join(errDir, ofn))
    mycmd = cmd + " >%s 2>%s" % (out_full_fn, err_full_fn)
    wd = os.path.join(workingDir, str(taski))
    mkdir_quiet(wd, 900)  # make the working directory
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
    for _ in xrange(args.num_retries + 1):
        if args.line_by_line:
            pipe = subprocess.Popen(mycmd, bufsize=-1, stdin=subprocess.PIPE, shell=True, cwd=wd)
            pipe.stdin.write(inp + '\n')
            pipe.stdin.close()
            ret = pipe.wait()
            if pipe.wait() == 0:
                return out_full_fn
            message('Non-zero return (%d) after closing pipe "%s"' % (ret, mycmd))
        else:
            ret = os.system(fullcmd)
            if ret == 0:
                return out_full_fn
            message('Non-zero return (%d) after executing "%s"' % (ret, fullcmd))
        message('Retrying in %d seconds...' % args.delay)
        time.sleep(args.delay)
        # Finished
    failQ.put([
        "Mapper %d of %d (pid %d) failed the maximum # of times %d" % (taski, taskn, os.getpid(), args.num_retries + 1),
        inp, err_full_fn, fullcmd])
    return out_full_fn


num_processes = min(num_processes, len(inps))

tasks, taski = [], 1
if args.line_by_line:
    for inp in inps:
        with open(inp) as fh:
            while True:
                ln = fh.readline().rstrip()
                if len(ln) == 0:
                    break
                if ln[0] == '#':
                    continue
                name = "%s:%d" % (inp, fh.tell())
                tasks.append((name, '\t'.join([name, ln]), taski))
                taski += 1
else:
    for inp in map(os.path.abspath, inps):
        tasks.append((inp, inp, taski))
        taski += 1
taskn = taski - 1

message('Piping %d task(s) to command "%s"' % (taskn, cmd))

pool = multiprocessing.Pool(num_processes)
outfns = []
r = pool.map_async(do_mapper, tasks, callback=outfns.extend)
while not r.ready():
    r.wait(1)

check_fail_queue()

if args.multiple_outputs:

    def do_split(tup, keep=args.keep_all):
        task, out_fn = tup
        splitdir = os.path.join(output, '_'.join([task, 'split']))
        check_dir(splitdir, 900)
        k2fh, k2fn = {}, {}
        with openex(out_fn) as fh:
            for ln in fh:
                ln = ln.rstrip()
                if len(ln) == 0: continue
                toks = string.split(ln, '\t')
                if toks[0] not in k2fh:
                    fn = os.path.join(splitdir, '_'.join([toks[0], task]))
                    k2fn[toks[0]] = fn
                    k2fh[toks[0]] = open(fn, 'w')
                k2fh[toks[0]].write('\t'.join(toks[1:]))
                k2fh[toks[0]].write('\n')
        for ofh in k2fh.itervalues():
            ofh.close()
        if not keep:
            os.remove(out_fn)
        return k2fn

    split_tups = []
    i = 0
    for outfn in outfns:
        split_tups.append((str(i), outfn))
        i += 1
    split_pool = multiprocessing.Pool(num_processes)
    k2fns = []
    r = split_pool.map_async(do_split, split_tups, callback=k2fns.extend)
    while not r.ready():
        r.wait(1)
    check_fail_queue()
    k2fn_list = defaultdict(list)
    for kToFn in k2fns:
        for k, v in kToFn.iteritems():
            k2fn_list[k].append(v)

    # Now join them back up
    for k, vl in k2fn_list.iteritems():
        k_out_dir = os.path.join(output, k)
        check_dir(k_out_dir, 900)
        for v in vl:
            fn = os.path.basename(v)
            shutil.copyfile(v, os.path.join(k_out_dir, fn))
            if not args.keep_all:
                os.remove(v)

    # Remove all the split directories
    if not args.keep_all:
        removed_already = set()
        for k, vl in k2fn_list.iteritems():
            for v in vl:
                vdir = os.path.dirname(v)
                if vdir not in removed_already:
                    shutil.rmtree(vdir)
                    removed_already.add(vdir)

if not args.keep_all:
    message('Removing intermediate directory "%s"\n' % intermediate)
    shutil.rmtree(intermediate)

message('SUCCESS\n')
