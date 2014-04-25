#!/usr/bin/perl

"""
reduce.py

Simple wrapper that mimics some of Hadoop's behavior during the Reduce step of a
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

parser = argparse.ArgumentParser(
    description='Simple wrapper that mimics some of Hadoop\'s behavior during the Map step of a MapReduce computation.')

parser.add_argument('--name', metavar='PATH', type=str, required=True, help='Name of this step')
parser.add_argument('--input', metavar='PATH', type=str, required=True, help='Files with input data')
parser.add_argument('--output', metavar='PATH', type=str, required=True, help='Files with output data')
parser.add_argument('--num-tasks', metavar='INT', type=int, required=True, help='Divide input into this many tasks.')
parser.add_argument('--bin-fields', metavar='INT', type=int, help='# fields to bin by.')
parser.add_argument('--sort-fields', metavar='INT', type=int, help='# fields to sort by.')
parser.add_argument('--external-sort', action='store_const', const=True, default=False,
                    help='Use external program to sort tuples.')
parser.add_argument('--sort-size', metavar='INT', type=int, default=1024 * 300, help='Memory cap on sorts.')
parser.add_argument('--messages', metavar='PATH', type=str, help='File to store stderr messages to')
parser.add_argument('--intermediate', metavar='PATH', type=str, help='Directory to store intermediate data in')
parser.add_argument('--write-every', type=int, default=200000000, help='Writes to disk only after this many characters have accumulated in memory')
parser.add_argument('--num-processes', metavar='INT', type=int,
                    help='Max # of simultaneous processes to run.  Default = # processors you have.')
parser.add_argument('--num-retries', metavar='INT', type=int, default=3,
                    help='# times to retry the mapper if something goes wrong')
parser.add_argument('--delay', metavar='INT', type=int, default=5, help='# seconds to wait before a retry')
parser.add_argument('--force', action='store_const', const=True, default=False, help='Profile the code')
parser.add_argument('--multiple-outputs', action='store_const', const=True, default=False,
                    help='Outputs go in subdirectories labeled with first output token')
parser.add_argument('--keep-all', action='store_const', const=True, default=False, help='Keep all intermediate results')
parser.add_argument('--profile', action='store_const', const=True, default=False, help='Profile the code')
parser.add_argument('--verbose', action='store_const', const=True, default=False,
                    help='Prints out extra debugging statements')

# Collect the reducer arguments
argv = sys.argv
reduce_argv = []
in_args = False
for i in xrange(1, len(sys.argv)):
    if in_args:
        reduce_argv.append(sys.argv[i])
    elif sys.argv[i] == '--':
        argv = sys.argv[:i]
        in_args = True

args = parser.parse_args(argv[1:])
assert args.num_tasks > 0

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

inps = args.input.split(',')
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
for inp in inps:
    message('  "%s"' % inp)
message('Output: "%s"' % output)
message('Intermediate: "%s"' % intermediate)
num_processes = args.num_processes or multiprocessing.cpu_count()
message('# parallel processes: %d' % num_processes)
message('Retries=%d, delay=%d seconds' % (args.num_retries, args.delay))
message('# tasks=%d' % args.num_tasks)
if args.num_tasks == 1:
    message('No binning and sorting')
else:
    message('Bin fields=%d, sort fields=%d' % (args.bin_fields, args.sort_fields))
message('Max sort memory footprint=%d' % args.sort_size)

options = []
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
err_dir = os.path.join(intermediate, 'reduce.err')
check_dir(err_dir, 300)

# Where reducer programs are actually run
working_dir = os.path.join(intermediate, 'reduce.wds')
check_dir(working_dir, 400)

# Where stderr output from sort processes go
sort_err_dir = os.path.join(intermediate, 'sort.err')
check_dir(sort_err_dir, 500)

# Where unsorted reduce-task inputs are stored temporarily
task_dir = os.path.join(intermediate, 'reduce.tasks')
check_dir(task_dir, 600)

# Where sorted reduce task inputs are stored temporarily
stask_dir = os.path.join(intermediate, 'reduce.stasks')
check_dir(stask_dir, 700)


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

fail_q = multiprocessing.Queue()


def check_fail_queue():
    if not fail_q.empty():
        while not fail_q.empty():
            err, inp, errfn, cmd = fail_q.get()
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

import random
message('Writing input records to %d tasks' % args.num_tasks)
task_names = [str(i) for i in xrange(args.num_tasks)]
task_lists = defaultdict(list)
tot = 0
char_count = 0
ofhs = {}
for inp in inps:
    with openex(inp) as fh:
        for ln in fh:
            if char_count > args.write_every:
                for task in task_lists:
                    task_lists[task].append('')
                    try:
                        ofhs[task].write('\n'.join(task_lists[task]))
                    except KeyError:
                        ofhs[task] = open(os.path.join(task_dir, task), 'w')
                        ofhs[task].write('\n'.join(task_lists[task]))
                char_count = 0
                task_lists = defaultdict(list)
            ln = ln.rstrip()
            ln_length = len(ln)
            if ln_length == 0:
                continue
            toks = string.split(ln, '\t')
            if toks[0] == 'DUMMY':
                continue
            tot += 1
            assert len(toks) >= args.bin_fields
            k = '\t'.join(toks[:args.bin_fields])
            random.seed(k)
            task = random.choice(task_names)
            task_lists[task].append(ln)
for task in task_lists:
    task_lists[task].append('')
    try:
        ofhs[task].write('\n'.join(task_lists[task]))
    except KeyError:
        ofhs[task] = open(os.path.join(task_dir, task), 'w')
        ofhs[task].write('\n'.join(task_lists[task]))
for fh in ofhs.itervalues():
    fh.close()
message('%d input records written' % tot)

'''Code below is memory-intensive if there are too many bins; seed
random number generator with key and choose random task instead, as above.
Note that this technically makes for more load imbalance among tasks, but 
it's less an issue in local mode, where Python processes will simply share 
resources as managed by the OS.'''
'''# Go through all input files in parallel, calculating bin sizes for all bins.
bin_count_queue = multiprocessing.Queue()


def bin_count_worker(fn):
    cnt = defaultdict(int)
    with openex(fn) as fh:
        for ln in fh:
            ln = ln.rstrip()
            if len(ln) == 0:
                continue
            toks = string.split(ln, '\t')
            if toks[0] == 'DUMMY':
                continue
            assert len(toks) >= args.bin_fields
            cnt['\t'.join(toks[:args.bin_fields])] += 1
    bin_count_queue.put(cnt)

bin_count = defaultdict(int)
count_pool = multiprocessing.Pool(num_processes)
r = count_pool.map_async(bin_count_worker, inps)
while not r.ready():
    r.wait(1)
check_fail_queue()

tot = 0
for _ in xrange(len(inps)):
    cnt = bin_count_queue.get()
    assert cnt is not None
    for k, v in cnt.iteritems():
        bin_count[k] += v
        tot += v

if tot > 0:
    message('Allocating %d input records to %d tasks' % (tot, args.num_tasks))

    # Allocate bins to tasks, always adding to task with least tuples so far
    task_names = ["task-%05d" % i for i in xrange(args.num_tasks)]
    taskq = [(0, task_names[i], []) for i in xrange(args.num_tasks)]
    key2task = {}
    for k, v in bin_count.iteritems():
        sz, nm, klist = heapq.heappop(taskq)
        sz += v
        klist.append(k)
        assert k not in key2task
        heapq.heappush(taskq, (sz, nm, klist))
        key2task[k] = nm

    message('Writing tasks')

    # Write out all the tasks to files within 'taskDir'
    ofhs = {}

    for inp in inps:
        with openex(inp) as fh:
            for ln in fh:
                ln = ln.rstrip()
                if len(ln) == 0:
                    continue
                toks = string.split(ln, '\t')
                if toks[0] == 'DUMMY':
                    continue
                assert len(toks) >= args.bin_fields
                k = '\t'.join(toks[:args.bin_fields])
                assert k in key2task, 'No such key: "%s"' % k
                task = key2task[k]
                if task not in ofhs:
                    ofhs[task] = open(os.path.join(task_dir, task), 'w')
                ofhs[task].write(ln)
                ofhs[task].write('\n')

    for fh in ofhs.itervalues():
        fh.close()'''

if tot > 0:

    ########################################
    # Stage 2. Sort and reduce each task
    ########################################

    message('Sorting tasks')
    need_sort = args.num_tasks > 1 or args.sort_fields > args.bin_fields
    if need_sort:
        def do_sort(task, external=True, keep=args.keep_all):
            assert external  # only know how to use external sort for now
            input_fn, sorted_fn = os.path.join(task_dir, task), os.path.join(stask_dir, task)
            sort_err_fn = os.path.join(sort_err_dir, task)
            cmd = 'sort -S %d -k1,%d %s >%s 2>%s' % (args.sort_size, args.sort_fields, input_fn, sorted_fn, sort_err_fn)
            el = os.system(cmd)
            if el != 0:
                msg = 'Sort command "%s" for sort task "%s" failed with exitlevel: %d' % (cmd, task, el)
                fail_q.put((msg, input_fn, sort_err_fn, cmd))
            elif not keep:
                os.remove(input_fn)
                os.remove(sort_err_fn)

        sortPool = multiprocessing.Pool(num_processes)
        r = sortPool.map_async(do_sort, ofhs.keys())
        while not r.ready():
            r.wait(1)
        check_fail_queue()

    reduce_cmd = ' '.join(reduce_argv)
    reduce_input_dir = stask_dir if need_sort else task_dir

    def do_reduce(task, keep=args.keep_all):
        sorted_fn = os.path.join(reduce_input_dir, task)
        sorted_fn = os.path.abspath(sorted_fn)
        if not os.path.exists(sorted_fn):
            raise RuntimeError('No such sorted task: "%s"' % sorted_fn)
        out_fn, err_fn = os.path.join(output, task), os.path.join(err_dir, task)
        out_fn, err_fn = os.path.abspath(out_fn), os.path.abspath(err_fn)
        cmd = 'cat %s | %s >%s 2>%s' % (sorted_fn, reduce_cmd, out_fn, err_fn)
        wd = os.path.join(working_dir, task)
        check_dir(wd, 800)
        message('Pid %d processing task %s' % (os.getpid(), task))
        # Simulate MapReduce task partition assignment
        new_env = os.environ.copy()
        new_env['mapred_task_partition'] = task
        pipe = subprocess.Popen(cmd, bufsize=-1, shell=True, cwd=wd, env=new_env)
        el = pipe.wait()
        if el != 0:
            msg = 'Reduce command "%s" for sort task "%s" failed with exitlevel: %d' % (cmd, task, el)
            fail_q.put((msg, sorted_fn, err_fn, cmd))
        elif not keep:
            os.remove(sorted_fn)
            os.remove(err_fn)
            shutil.rmtree(wd)
        return out_fn

    message('Piping %d task(s) to command "%s"' % (len(ofhs.keys()), reduce_cmd))
    reduce_pool = multiprocessing.Pool(num_processes)
    outfns = []
    r = reduce_pool.map_async(do_reduce, ofhs.keys(), callback=outfns.extend)
    while not r.ready():
        r.wait(1)
    check_fail_queue()

    if args.multiple_outputs:

        def do_split(task, keep=args.keep_all):
            out_fn = os.path.join(output, task)
            split_dir = os.path.join(output, '_'.join([task, 'split']))
            check_dir(split_dir, 900)
            k2fh, k2fn = {}, {}
            with openex(out_fn) as fh:
                for ln in fh:
                    ln = ln.rstrip()
                    if len(ln) == 0: continue
                    toks = string.split(ln, '\t')
                    if toks[0] not in k2fh:
                        fn = os.path.join(split_dir, '_'.join([toks[0], task]))
                        k2fn[toks[0]] = fn
                        k2fh[toks[0]] = open(fn, 'w')
                    k2fh[toks[0]].write('\t'.join(toks[1:]))
                    k2fh[toks[0]].write('\n')
            for ofh in k2fh.itervalues():
                ofh.close()
            if not keep:
                os.remove(out_fn)
            return k2fn

        split_pool = multiprocessing.Pool(num_processes)
        k2fns = []
        r = split_pool.map_async(do_split, ofhs.keys(), callback=k2fns.extend)
        while not r.ready():
            r.wait(1)
        check_fail_queue()
        k2fn_list = defaultdict(list)
        for k2fn in k2fns:
            for k, v in k2fn.iteritems():
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
else:
    message("0 input records!  Skipping binning, sorting, reducing, etc...")

if not args.keep_all:
    message('Removing intermediate directory "%s"\n' % intermediate)
    shutil.rmtree(intermediate)

message('SUCCESS\n')
