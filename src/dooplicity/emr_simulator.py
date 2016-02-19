#!/usr/bin/env python
"""
emr_simulator.py
Part of Dooplicity framework

Runs JSON-encoded Hadoop Streaming job flow. FUNCTIONALITY IS IDIOSYNCRATIC;
it is currently confined to those features used by Rail. Format of input JSON
mirrors that of StepConfig list from JSON sent to EMR via RunJobsFlow. Any
files input to a mapper can be gzip'd, but inputs to a reducer currently cannot
be.

In --ipy mode, the script uses IPython to run tasks on different engines
mediated by a controller. IPython controller and engines must be started before
this script is invoked.

All paths in input JSON should be absolute.

Licensed under the MIT License:

Copyright (c) 2014 Abhi Nellore and Ben Langmead.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import argparse
import sys
from collections import defaultdict, OrderedDict, deque
import time
import json
import interface as dp_iface
import gc
import shutil
import signal
import itertools
import socket
import subprocess
import glob
import hashlib
import tempfile
import shutil
import os
import contextlib
from tools import make_temp_dir, make_temp_dir_and_register_cleanup
from ansibles import Url
import site
import string

def add_args(parser):
    """ Adds args relevant to EMR simulator.

        parser: object of type parser.ArgumentParser

        No return value.
    """
    parser.add_argument(
            '-m', '--memcap', type=int, required=False, default=(1024*300),
            help=('Maximum amount of memory (in bytes) to use per UNIX sort '
                  'instance.')
        )
    parser.add_argument(
            '-p', '--num-processes', type=int, required=False, default=1,
            help='Number of subprocesses to open at once.'
        )
    parser.add_argument(
            '-t', '--max-attempts', type=int, required=False, default=4,
            help=('Maximum number of times to attempt a task.')
        )
    parser.add_argument(
            '-s', '--separator', type=str, required=False, default='\t',
            help='Separator between successive fields in inputs and '
                 'intermediates.'
        )
    parser.add_argument(
            '-k', '--keep-intermediates', action='store_const', const=True,
            default=False,
            help='Keeps all intermediate output.'
        )
    parser.add_argument(
            '--keep-last-output', action='store_const', const=True,
            default=False,
            help='If --keep-intermediates is False, keeps outputs that are ' \
                 'unused as inputs by steps.'
        )
    parser.add_argument('--gzip-outputs', action='store_const',
            const=True, default=False,
            help='Compress step output files with gzip.'
        )
    parser.add_argument('--gzip-level', type=int, required=False,
            default=3,
            help='Level of gzip compression to use, if applicable.'
        )
    parser.add_argument('--ipy', action='store_const', const=True,
            default=False,
            help=('Uses IPython controller and engines to execute tasks; this '
                  'permits running a MapReduce job flow on a wide array of '
                  'cluster setups. Ignores --num-processes in favor of the '
                  'number of available engines.')
        )
    parser.add_argument('--ipcontroller-json', type=str, required=False,
            default=None,
            help=('Path to ipcontroller-client.json file; relevant only if '
                  '--ipy is invoked. See IPython documentation for '
                  'more information. If left unspecified, IPython\'s default '
                  'path is used.')
        )
    parser.add_argument('--ipy-profile', type=str, required=False,
            default=None,
            help=('Connects to this IPython profile; relevant only if --ipy '
                  'is invoked and takes precedence over --ipcontroller-json.')
        )
    parser.add_argument('--scratch', type=str, required=False,
            default=None,
            help=('Where to write any intermediate output before copying to '
                  'consolidated intermediate directory. This is typically '
                  'a directory local to a given node. None means write '
                  'directory to consolidated intermediate directory. The '
                  'string \"-\" means write to a temporary directory securely '
                  'created by Python.')
        )
    parser.add_argument('--direct-write', action='store_const',
            const=True, default=False,
            help=('Always write intermediate files directly to consolidated '
                  'intermediate directory, even if --scratch is specified.')
        )
    parser.add_argument('--common', type=str, required=False,
            default=None,
            help=('Location of a writable directory accessible across all '
                  'nodes; this is where some temporary files may be stored '
                  'and is not important unless running in --ipy mode; if '
                  'left unspecified, defaults to Python temporary directory'))
    parser.add_argument('--sort', type=str, required=False,
            default='sort',
            help=('Path to sort executable. Add arguments as necessary, '
                  'e.g. for specifying a directory for storing sort\'s '
                  'temporary files.'))

def init_worker():
    """ Prevents KeyboardInterrupt from reaching a pool's workers.

        Exiting gracefully after KeyboardInterrupt or SystemExit is a
        challenge. The solution implemented here is by John Reese and is from
        http://noswap.com/blog/python-multiprocessing-keyboardinterrupt .

        No return value.
    """
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def yopen(gzipped, *args):
    """ Passes args on to the appropriate opener, gzip or regular.

        A dooplicity.tools.xopen that uses the gzip module, which is
        unsafe for writing. See xopen's docstring in dooplicity.tools for
        more information.

        gzipped: True iff gzip.open() should be used to open rather than
            open(); False iff open() should be used; None if input should be
            read and guessed
        *args: unnamed arguments to pass

        Return value: file object
    """
    import gzip
    if gzipped is None:
        with open(args[0], 'rb') as binary_input_stream:
            # Check for magic number
            if binary_input_stream.read(2) == '\x1f\x8b':
                gzipped = True
            else:
                gzipped = False
    if gzipped:
        return gzip.open(*args)
    return open(*args)

def parsed_keys(partition_options, key_fields):
    """ Parses UNIX sort options to figure out what to partition on.

        Returned is a function that takes a line as input and returns a tuple
        of elements from the line to partition on OR False if the args are
        invalid.

        partition_options: UNIX sort options like -k1,1 -k3 -k3,4r -k 4 -k 5,3
        key_fields: number of fields from line to consider key

        Return value: see above
    """
    try:
        # Make list of tuples of start, end indexes
        parsed_args = [
                tuple([int(el) - 1 
                        for el in arg.strip().strip('nr').split(',')])
                        for arg in partition_options.split('-k')
                        if arg.strip() and len(arg.split(',')) <= 2
            ]
    except Exception:
        # args are invalid
        return False
    else:
        exec (
"""def partitioned_key(line, separator):
    key = line.strip().split(separator)[:{key_fields}]
    return {return_value}
""".format(key_fields=key_fields,
            return_value='+'.join(['key[{}:{}]'.format(
                                arg[0], arg[1] + 1 if len(arg) == 2 else ''
                            ) for arg in parsed_args]))
        )
        return partitioned_key

def presorted_tasks(input_files, process_id, sort_options, output_dir,
                    key_fields, separator, partition_options, task_count,
                    memcap, gzip=False, gzip_level=3, scratch=None,
                    direct_write=False, sort='sort', mod_partition=False,
                    max_attempts=4):
    """ Partitions input data into tasks and presorts them.

        Files in output directory are in the format x.y, where x is a task
        number on the interval [0, number of tasks - 1], and y is a process
        ID that identifies which process created the file. y is unimportant;
        the glob x.* should be catted to the reducer.

        Formula for computing task assignment: 
            int(hashlib.md5(key).hexdigest(), 16) % (task_count)

        input_files: list of files on which to operate.
        process_id: unique identifier for current process.
        sort_options: options to use when presorting.
        output_dir: directory in which to write output files.
        key_fields: number of fields from a line to consider the key.
        separator: separator between successive fields from line.
        partition_options: sort-like options to use when partitioning.
        streaming_command: streaming command to run.
        task_count: number of tasks in which to partition input.
        memcap: maximum percent of memory to use per UNIX sort instance.
        gzip: True iff all files written should be gzipped; else False.
        gzip_level: Level of gzip compression to use, if applicable.
        scratch: where to write output before copying to output_dir. If "-"
            string, writes to temporary directory; if None, writes directly
            to output directory.
        direct_write: write intermediate files directly to final destination,
            no matter what scratch is.
        sort: path to sort executable
        mod_partition: if True, task is assigned according to formula
            (product of fields) % task_count
        max_attempts: maximum number of times to attempt partitioning input.
            MUST BE FINAL ARG to be compatible with 
            execute_balanced_job_with_retries().

        Return value: None if no errors encountered; otherwise error string.
    """
    try:
        from operator import mul
        task_streams = {}
        if scratch is not None:
            scratch = os.path.expanduser(os.path.expandvars(scratch))
        if gzip:
            task_stream_processes = {}
        if direct_write:
            final_output_dir = output_dir
        elif scratch == '-':
            # Write to temporary directory
            final_output_dir = output_dir
            try:
                output_dir = tempfile.mkdtemp()
            except OSError as e:
                return ('Problem encountered creating temporary '
                        'scratch subdirectory: %s' % e)
        elif scratch:
            # Write to temporary directory in special location
            final_output_dir = output_dir
            try:
                os.makedirs(scratch)
            except OSError as e:
                if os.path.isfile(scratch):
                    return ('Scratch directory %s is a file.' % scratch)
            except IOError as e:
                return ('Scratch directory %s is not writable: %s' % (scratch,
                                                                        e))
            try:
                output_dir = tempfile.mkdtemp(dir=scratch)
            except OSError as e:
                return ('Problem encountered creating temporary '
                        'scratch subdirectory of %s: %s' % (scratch, e))
        else:
            final_output_dir = output_dir
        output_dir = os.path.expandvars(output_dir)
        final_output_dir = os.path.expandvars(final_output_dir)
        partitioned_key = parsed_keys(partition_options, separator)
        if not partitioned_key:
            # Invalid partition options
            return ('Partition options "%s" are invalid.' % partition_options)
        for input_file in input_files:
            with yopen(None, input_file) as input_stream:
                for line in input_stream:
                    key = partitioned_key(line, separator)
                    if mod_partition and len(key) <= 1:
                        try:
                            task = abs(int(key[0])) % task_count
                        except (IndexError, ValueError):
                            # Null key or some field doesn't work with this
                            task = int(
                                hashlib.md5(separator.join(key)).hexdigest(),
                                16
                            ) % task_count
                    else:
                        task = int(
                            hashlib.md5(separator.join(key)).hexdigest(), 16
                        ) % task_count
                    try:
                        task_streams[task].write(line)
                    except KeyError:
                        # Task file doesn't exist yet; create it
                        if gzip:
                            task_file = os.path.join(output_dir, str(task) +
                                                        '.' + str(process_id)
                                                        + '.unsorted.gz')
                            task_stream_processes[task] = subprocess.Popen(
                                    'gzip -%d >%s' % 
                                    (gzip_level, task_file),
                                    shell=True, bufsize=-1,
                                    executable='/bin/bash',
                                    stdin=subprocess.PIPE
                                )
                            task_streams[task] \
                                = task_stream_processes[task].stdin
                        else:
                            task_file = os.path.join(output_dir, str(task) +
                                                        '.' + str(process_id)
                                                        + '.unsorted')
                            task_streams[task] = open(task_file, 'w')
                        task_streams[task].write(line)
        for task in task_streams:
            task_streams[task].close()
        if gzip:
            for task in task_stream_processes:
                task_stream_processes[task].wait()
        # Presort task files
        if gzip:
            for unsorted_file in glob.glob(os.path.join(
                                                    output_dir,
                                                    '*.%s.unsorted.gz'
                                                    % process_id
                                                )):
                sort_command = (('set -eo pipefail; gzip -cd %s | '
                                 'LC_ALL=C %s -S %d %s -t$\'%s\' | '
                                 'gzip -c -%d >%s')
                                    % (unsorted_file, sort, memcap,
                                        sort_options,
                                        separator.encode('string_escape'),
                                        gzip_level,
                                        unsorted_file[:-12] + '.gz'))
                try:
                    subprocess.check_output(sort_command,
                                            shell=True,
                                            executable='/bin/bash',
                                            bufsize=-1,
                                            stderr=subprocess.STDOUT)
                except subprocess.CalledProcessError as e:
                    return (('Error "%s" encountered sorting file %s; exit '
                             'code was %d; command invoked was "%s".') %
                                (e.output.strip(),
                                    unsorted_file, e.returncode,
                                    sort_command))
                finally:
                    os.remove(unsorted_file)
        else:
            for unsorted_file in glob.glob(os.path.join(
                                                    output_dir,
                                                    '*.%s.unsorted'
                                                    % process_id
                                                )):
                sort_command = 'LC_ALL=C %s -S %d %s -t$\'%s\' %s >%s' % (
                                                            sort, memcap,
                                                            sort_options,
                                                            separator.encode(
                                                                'string_escape'
                                                            ),
                                                            unsorted_file,
                                                            unsorted_file[:-9]
                                                        )
                try:
                    subprocess.check_output(sort_command,
                                              shell=True,
                                              executable='/bin/bash',
                                              bufsize=-1,
                                              stderr=subprocess.STDOUT)
                except subprocess.CalledProcessError as e:
                    return (('Error "%s" encountered sorting file %s; exit '
                             'code was %d; command invoked was "%s".') %
                                (e.output.strip(),
                                    unsorted_file, e.returncode,
                                    sort_command))
                finally:
                    os.remove(unsorted_file)
        return None
    except Exception:
        # Uncaught miscellaneous exception
        from traceback import format_exc
        return ('Error\n\n%s\nencountered partitioning input files '
                '[%s] into tasks.'
                        % (format_exc(),
                            (('%s, '* (len(input_files) - 1) 
                                       + '%s') % tuple(input_files))))
    finally:
        if 'final_output_dir' in locals() and final_output_dir != output_dir:
            # Copy all output files to final destination and kill temp dir
            for root, dirnames, filenames in os.walk(output_dir):
                if not filenames: continue
                destination = os.path.join(
                                    final_output_dir,
                                    os.path.relpath(root, output_dir)
                                )
                try:
                    os.makedirs(destination)
                except OSError:
                    # Directory already exists
                    pass
                for filename in filenames:
                    shutil.copy(
                            os.path.join(root, filename),
                            os.path.join(destination, filename)
                        )
            shutil.rmtree(output_dir)

def step_runner_with_error_return(streaming_command, input_glob, output_dir,
                                  err_dir, task_id, multiple_outputs,
                                  separator, sort_options, memcap,
                                  gzip=False, gzip_level=3, scratch=None,
                                  direct_write=False, sort='sort',
                                  dir_to_path=None, attempt_number=None):
    """ Runs a streaming command on a task, segregating multiple outputs. 

        streaming_command: streaming command to run.
        input_glob: input files on which to run streaming command
            specified with wildcard; files in a directory are
            considered, while subdirectories are neglected. Either every
            file in input_glob should be gzip'd, or none should be.
            Gzip'd input is accommodated only for a map step!
        output_dir: directory in which to write output.
        err_dir: directory in which to write errors
        task_id: unique numerical identifer for task. Used to set
            environment variable mapred_task_partition and determine
            output filename.
        multiple_outputs: True if output should be divided by key before
            first instance of separator, described below.
        separator: character separating successive fields in a line from
            input_file.
        sort_options: None if no sort should be performed on input_glob;
            otherwise, performs merge sort with unix sort -m and the
            specified string of command-line parameters. EACH INPUT FILE
            SHOULD BE PRESORTED.
        memcap: maximum percent of memory to use per UNIX sort instance.
        gzip: True iff all files written should be gzipped; else False.
        gzip_level: Level of gzip compression to use, if applicable.
        scratch: where to write output before copying to output_dir. If "-"
            string, writes to temporary directory; if None, writes directly
            to output directory.
        direct_write: write intermediate files directly to final destination,
            no matter what scratch is.
        sort: path to sort executable.
        dir_to_path: path to add to PATH.
        attempt_number: attempt number of current task or None if no retries.
            MUST BE FINAL ARG to be compatible with 
            execute_balanced_job_with_retries().

        Return value: None iff step runs successfully; otherwise error message.
    """
    command_to_run = None
    try:
        if direct_write:
            final_output_dir = output_dir
        elif scratch == '-':
            # Write to temporary directory
            final_output_dir = output_dir
            try:
                output_dir = tempfile.mkdtemp()
            except OSError as e:
                return ('Problem encountered creating temporary '
                        'scratch subdirectory: %s' % e)
        elif scratch:
            scratch = os.path.expanduser(os.path.expandvars(scratch))
            # Write to temporary directory in special location
            final_output_dir = output_dir
            try:
                os.makedirs(scratch)
            except OSError as e:
                if os.path.isfile(scratch):
                    return ('Scratch directory %s is a file.' % scratch)
            except IOError as e:
                return ('Scratch directory %s is not writable: %s' % (scratch,
                                                                        e))
            try:
                output_dir = tempfile.mkdtemp(dir=scratch)
            except OSError as e:
                return ('Problem encountered creating temporary '
                        'scratch subdirectory of %s: %s' % (scratch, e))
        else:
            final_output_dir = output_dir
        output_dir = os.path.expandvars(output_dir)
        final_output_dir = os.path.expandvars(final_output_dir)
        input_files = [input_file for input_file in glob.glob(input_glob)
                            if os.path.isfile(input_file)]
        if not input_files:
            # No input!
            return None
        if sort_options is None:
            # Mapper. Check if first input file is gzip'd
            with open(input_files[0], 'rb') as binary_input_stream:
                if binary_input_stream.read(2) == '\x1f\x8b':
                    # Magic number of gzip'd file found
                    prefix = 'gzip -cd %s' % input_glob
                else:
                    prefix = 'cat %s' % input_glob
        else:
            # Reducer. Merge sort the input glob.
            if gzip:
                # Use process substitution
                prefix = '(LC_ALL=C %s -S %d %s -t$\'%s\' -m %s' % (
                        sort, memcap, sort_options,
                        separator.encode('string_escape'),
                        ' '.join(['<(gzip -cd %s)' % input_file
                                    for input_file in input_files]) + ')'
                    )
            else:
                # Reducer. Merge sort the input glob.
                prefix = 'LC_ALL=C %s -S %d %s -t$\'%s\' -m %s' % (sort,
                                            memcap,
                                            sort_options,
                                            separator.encode('string_escape'),
                                            input_glob)
        err_file = os.path.abspath(os.path.join(err_dir, (
                                            ('%d.log' % task_id)
                                                if attempt_number is None
                                                else ('%d.%d.log'
                                                    % (task_id, attempt_number)
                                                )
                                            )
                                        ))
        new_env = os.environ.copy()
        new_env['mapreduce_task_partition'] \
            = new_env['mapred_task_partition'] = str(task_id)
        if multiple_outputs:
            # Must grab each line of output and separate by directory
            command_to_run \
                = prefix + ' | ' + streaming_command + (' 2>%s' % err_file)
            # Need bash or zsh for process substitution
            multiple_output_process = subprocess.Popen(
                    ' '.join([('set -eo pipefail; cd %s;' % dir_to_path)
                                if dir_to_path is not None
                                else 'set -eo pipefail;',
                              command_to_run]),
                    shell=True,
                    stdout=subprocess.PIPE,
                    stderr=open(os.devnull, 'w'),
                    env=new_env,
                    bufsize=-1,
                    executable='/bin/bash'
                )
            task_file_streams = {}
            if gzip:
                task_file_stream_processes = {}
            for line in multiple_output_process.stdout:
                key, _, line_to_write = line.partition(separator)
                try:
                    task_file_streams[key].write(line_to_write)
                except KeyError:
                    '''Must create new file, but another process could have
                    created the output directory.'''
                    key_dir = os.path.join(output_dir, key)
                    try:
                        os.makedirs(key_dir)
                    except OSError:
                        if not os.path.exists(key_dir):
                            return (('Streaming command "%s" failed: problem '
                                     'encountered creating output '
                                     'directory %s.') % (command_to_run,
                                                          key_dir))
                    if gzip:
                        task_file_stream_processes[key] = subprocess.Popen(
                                'gzip -%d >%s' % 
                                (gzip_level,
                                 os.path.join(key_dir, str(task_id) + '.gz')),
                                shell=True, bufsize=-1,
                                executable='/bin/bash',
                                stdin=subprocess.PIPE
                            )
                        task_file_streams[key] \
                            = task_file_stream_processes[key].stdin
                    else:
                        task_file_streams[key] = open(
                                os.path.join(key_dir, str(task_id)), 'w'
                            )
                    task_file_streams[key].write(line_to_write)
            multiple_output_process_return = multiple_output_process.wait()
            if multiple_output_process_return != 0:
                return (('Streaming command "%s" failed; exit level was %d.')
                         % (command_to_run, multiple_output_process_return))
            for key in task_file_streams:
                task_file_streams[key].close()
        else:
            if gzip:
                out_file = os.path.abspath(
                                os.path.join(output_dir, str(task_id) + '.gz')
                            )
                command_to_run \
                    = prefix + ' | ' + streaming_command + (
                            ' 2>%s | gzip -%d >%s'
                                % (err_file,
                                    gzip_level,
                                    out_file)
                        )
            else:
                out_file = os.path.abspath(
                                os.path.join(output_dir, str(task_id))
                            )
                command_to_run \
                    = prefix + ' | ' + streaming_command + (' >%s 2>%s'
                                                             % (out_file,
                                                                err_file))
            try:
                # Need bash or zsh for process substitution
                subprocess.check_output(' '.join([('set -eo pipefail; cd %s;'
                                                    % dir_to_path)
                                                    if dir_to_path is not None
                                                    else 'set -eo pipefail;',
                                                  command_to_run]),
                                            shell=True,
                                            env=new_env,
                                            bufsize=-1,
                                            stderr=subprocess.STDOUT,
                                            executable='/bin/bash')
            except subprocess.CalledProcessError as e:
                return (('Streaming command "%s" failed; exit level was %d.')
                         % (command_to_run, e.returncode))
        return None
    except Exception as e:
        # Uncaught miscellaneous exception
        from traceback import format_exc
        return ('Error\n\n%s\nencountered executing task on input %s.'
                % (format_exc(), input_glob))
    finally:
        if 'task_file_stream_processes' in locals():
            for key in task_file_streams:
                task_file_streams[key].close()
            for key in task_file_stream_processes:
                task_file_stream_processes[key].wait()
        if 'final_output_dir' in locals() and final_output_dir != output_dir:
            # Copy all output files to final destination and kill temp dir
            for root, dirnames, filenames in os.walk(output_dir):
                if not filenames: continue
                destination = os.path.join(
                                    final_output_dir,
                                    os.path.relpath(root, output_dir)
                                )
                try:
                    os.makedirs(destination)
                except OSError:
                    # Directory already exists
                    pass
                for filename in filenames:
                    shutil.copy(
                            os.path.join(root, filename),
                            os.path.join(destination, filename)
                        )
            shutil.rmtree(output_dir)

def run_simulation(branding, json_config, force, memcap, num_processes,
                    separator, keep_intermediates, keep_last_output,
                    log, gzip=False, gzip_level=3, ipy=False,
                    ipcontroller_json=None, ipy_profile=None, scratch=None,
                    common=None, sort='sort', max_attempts=4,
                    direct_write=False):
    """ Runs Hadoop Streaming simulation.

        FUNCTIONALITY IS IDIOSYNCRATIC; it is currently confined to those
        features used by Rail. Format of input JSON mirrors that used by
        elastic-mapreduce-ruby. Any files input to a mapper can be gzip'd,
        but inputs to a reducer currently cannot be.

        branding: text file with branding to print to screen before running
            job. This is where the name of a software package or ASCII art 
            can go.
        json_config: JSON configuration file. Google Getting Started with
            Amazon Elastic MapReduce for formatting information.
        force: True iff all existing directories should be erased when
            writing intermediates.
        memcap: maximum fraction of memory to use across UNIX sort instances.
        num_processes: number of subprocesses to open at once; applicable only 
            when not in ipy mode
        separator: separator between successive fields in inputs and 
            intermediates.
        keep_intermediates: keeps all intermediate output.
        keep_last_output: keeps outputs that are unused as inputs by steps.
        log: name of file in which to store messages written to stderr
        gzip: True iff all files written should be gzipped; else False.
        gzip_level: level of gzip compression to use, if applicable.
        ipy: use iPython engines to run tasks.
        ipcontroller_json: path to ipcontroller-client.json; relevant only if 
            ipy is True. If None, uses IPython's default location.
        ipy_profile: name of IPython cluster configuration profile to use; None
            if profile is not specified. In this case, ipcontroller_json takes
            precedence.
        scratch: scratch directory, typically local. Files are written here by
            tasks before copying to final destination. If None, files are
            written directly to final destination. If '-', files are written
            to safely created temporary directory.
        common: path to directory accessible across nodes in --ipy mode
        sort: sort executable including command-line arguments
        max_attempts: maximum number of times to attempt a task in ipy mode.
        direct_write: always writes intermediate files directly to final
            destination, even when scratch is specified

        No return value.
    """
    global failed
    import shutil
    import os
    import tempfile
    import glob
    if log is not None:
        try:
            os.makedirs(os.path.dirname(log))
        except OSError:
            pass
        try:
            log_stream = open(log, 'a')
        except Exception as e:
            log_stream = None
    else:
        log_stream = None
    iface = dp_iface.DooplicityInterface(branding=branding,
                                         log_stream=log_stream)
    failed = False
    try:
        # Using IPython?
        if ipy:
            try:
                from IPython.parallel import Client
            except ImportError:
                iface.fail('IPython is required to run Dooplicity\'s EMR '
                           'simulator in --ipy mode. Visit ipython.org to '
                           'download it, or simply download the Anaconda '
                           'distribution of Python at '
                           'https://store.continuum.io/cshop/anaconda/; it\'s '
                           'easy to install and comes with IPython and '
                           'several other useful packages.')
                failed = True
                raise
            if ipy_profile:
                try:
                    pool = Client(profile=ipy_profile)
                except ValueError:
                    iface.fail('Cluster configuration profile "%s" was not '
                               'found.' % ipy_profile)
                    failed = True
                    raise
            elif ipcontroller_json:
                try:
                    pool = Client(ipcontroller_json)
                except IOError:
                    iface.fail(
                            'Cannot find connection information JSON file %s.'
                            % ipcontroller_json
                        )
                    failed = True
                    raise
            else:
                try:
                    pool = Client()
                except IOError:
                    iface.fail(
                            'Cannot find ipcontroller-client.json. Ensure '
                            'that IPython controller and engines are running.'
                            ' If controller is running on a remote machine, '
                            'copy the ipcontroller-client.json file from there '
                            'to a local directory; then rerun this script '
                            'specifying the local path to '
                            'ipcontroller-client.json with the '
                            '--ipcontroller-json command-line parameter.'
                        )
                    failed = True
                    raise
            if not pool.ids:
                iface.fail(
                        'An IPython controller is running, but no engines are '
                        'connected to it.'
                    )
                failed = True
                raise RuntimeError
            # Use all engines
            num_processes = len(pool)
            all_engines = set(pool.ids)
            from tools import apply_async_with_errors
            direct_view = pool[:]
            # Use Dill to permit general serializing
            try:
                import dill
            except ImportError:
                raise RuntimeError(
                        'Dooplicity requires Dill. Install it by running '
                        '"pip install dill", or see the StackOverflow '
                        'question http://stackoverflow.com/questions/23576969/'
                        'how-to-install-dill-in-ipython for other leads.'
                    )
            else:
                direct_view.use_dill()
            iface.status('Loading dependencies on IPython engines...')
            with direct_view.sync_imports(quiet=True):
                import subprocess
                import glob
                import hashlib
                import tempfile
                import shutil
                import os
            direct_view.push(dict(
                    yopen=yopen,
                    step_runner_with_error_return=\
                        step_runner_with_error_return,
                    presorted_tasks=presorted_tasks,
                    parsed_keys=parsed_keys
                ))
            iface.step('Loaded dependencies on IPython engines.')
            # Get host-to-engine and engine pids relations
            current_hostname = socket.gethostname()
            host_map = apply_async_with_errors(
                                    pool, all_engines, socket.gethostname,
                                    dict_format=True
                                )
            engine_map = defaultdict(list)
            for engine in host_map:
                engine_map[host_map[engine]].append(engine)
            pid_map = apply_async_with_errors(
                                    pool, all_engines, os.getpid,
                                    dict_format=True
                                )
            def interrupt_engines(pool, iface):
                """ Interrupts IPython engines spanned by view

                    Taken from:
                    http://mail.scipy.org/pipermail/ipython-dev/
                    2014-March/013426.html

                    pool: IPython Client object
                    iface: instance of DooplicityInterface

                    No return value.
                """
                iface.status('Interrupting IPython engines...')
                for engine_id in pool.ids:
                    host = host_map[engine_id]
                    kill_command = (
                          'CPIDS=$(pgrep -P {}); echo $CPIDS;'
                          '(sleep 33 && kill -9 $CPIDS &); '
                          'kill -9 $CPIDS'
                        ).format(pid_map[engine_id])
                    if host == socket.gethostname():
                        pass
                        # local
                        #subprocess.Popen(kill_command,
                        #        bufsize=-1, shell=True
                        #    )
                    else:
                        #subprocess.Popen(
                        #    ('ssh -oStrictHostKeyChecking=no '
                        #     '-oBatchMode=yes {} \'{}\'').format(
                        #        host, kill_command
                        #    ), bufsize=-1, shell=True
                        #)
                        pass
            import random
            def execute_balanced_job_with_retries(pool, iface,
                task_function, task_function_args,
                status_message='Tasks completed',
                finish_message='Completed tasks.', max_attempts=4):
                """ Executes parallel job over IPython engines with retries.

                    Tasks are assigned to free engines as they become
                    available. If a task fails on one engine, it is retried on
                    another engine. If a task has been tried on all engines but
                    fails before max_attempts is exceeded, the step is failed.

                    pool: IPython Client object; all engines it spans are used
                    iface: DooplicityInterface object for spewing log messages
                        to console
                    task_function: name if function to execute
                    task_function_args: iterable of lists, each of whose
                        items are task_function's arguments, WITH THE EXCEPTION
                        OF A SINGLE KEYWORD ARGUMENT "attempt_count". This
                        argument must be the final keyword argument of the
                        function but _excluded_ from the arguments in any item
                        of task_function_args.
                    status_message: status message about tasks completed
                    finish_message: message to output when all tasks are
                        completed
                    max_attempts: max number of times to attempt any given
                        task

                    No return value.
                """
                global failed
                random.seed(pool.ids[-1])
                used_engines, free_engines = set(), set(pool.ids)
                completed_tasks = 0
                tasks_to_assign = deque([
                        [task_function_arg, i, []] for i, task_function_arg
                        in enumerate(task_function_args)
                    ])
                task_count = len(tasks_to_assign)
                assigned_tasks, asyncresults = {}, {}
                max_task_fails = 0
                iface.status(('    %s: '
                              '%d/%d | \\max_i (task_i fails): %d/%d')
                                % (status_message, completed_tasks,
                                    task_count, max_task_fails,
                                    max_attempts - 1))
                while completed_tasks < task_count:
                    if tasks_to_assign:
                        task_to_assign = tasks_to_assign.popleft()
                        forbidden_engines = set(task_to_assign[2])
                        if len(forbidden_engines) >= 2:
                            # After two fails, do not allow reused nodes
                            for forbidden_engine in task_to_assign[2]:
                                forbidden_engines.update(
                                    engine_map[host_map[forbidden_engine]]
                                )
                        if all_engines <= forbidden_engines:
                            iface.fail(('No more running IPython engines '
                                        'and/or nodes on which function-arg '
                                        'combo (%s, %s) has not failed '
                                        'attempt to execute. Check the '
                                        'IPython cluster\'s integrity and '
                                        'resource availability.')
                                         % (task_function, task_to_assign[0]),
                                         steps=(job_flow[step_number:]
                                            if step_number != 0 else None))
                            failed = True
                            raise RuntimeError
                        try:
                            assigned_engine = random.choice(
                                    list(free_engines - forbidden_engines)
                                )
                        except IndexError:
                            # No engine to assign yet; add back to queue
                            tasks_to_assign.append(task_to_assign)
                        else:
                            asyncresults[task_to_assign[1]] = (
                                pool[assigned_engine].apply_async(
                                    task_function,
                                    *(task_to_assign[0] +
                                      [len(task_to_assign[2])])
                                )
                            )
                            assigned_tasks[task_to_assign[1]] = [
                                    task_to_assign[0], task_to_assign[1],
                                    task_to_assign[2] + [assigned_engine]
                                ]
                            used_engines.add(assigned_engine)
                            free_engines.remove(assigned_engine)
                    asyncresults_to_remove = []
                    for task in asyncresults:
                        if asyncresults[task].ready():
                            return_value = asyncresults[task].get()
                            if return_value is not None:
                                if max_attempts > len(assigned_tasks[task][2]):
                                    # Add to queue for reattempt
                                    tasks_to_assign.append(
                                            assigned_tasks[task]
                                        )
                                    max_task_fails = max(
                                            len(assigned_tasks[task][2]),
                                            max_task_fails
                                        )
                                    asyncresults_to_remove.append(task)
                                else:
                                    # Bail if max_attempts is saturated
                                    iface.fail(return_value,
                                    steps=(job_flow[step_number:]
                                            if step_number != 0 else None))
                                    failed = True
                                    raise RuntimeError
                            else:
                                # Success
                                completed_tasks += 1
                                asyncresults_to_remove.append(task)
                            iface.status(('    %s: '
                                          '%d/%d | '
                                          '\\max_i (task_i fails): '
                                          '%d/%d')
                                % (status_message, completed_tasks,
                                    task_count, max_task_fails,
                                    max_attempts - 1))
                            assert assigned_tasks[task][-1][-1] == \
                                asyncresults[task].engine_id
                            # Free engine
                            used_engines.remove(
                                    assigned_tasks[task][-1][-1]
                                )
                            free_engines.add(assigned_tasks[task][-1][-1])
                    for task in asyncresults_to_remove:
                        del asyncresults[task]
                        del assigned_tasks[task]
                    time.sleep(0.1)
                assert not used_engines
                iface.step(finish_message)
            @contextlib.contextmanager
            def cache(pool=None, file_or_archive=None, archive=True):
                """ Places X.[tar.gz/tgz]#Y in dir Y, unpacked if archive

                    pool: IPython Client object; all engines it spans are used
                    archive: file in format X.tar.gz#Y; None if nothing should
                        be done

                    Yields before deleting all files.
                """
                global failed
                if file_or_archive is None:
                    '''So with statements can be used all the time, even
                    when there's nothing to be archived'''
                    assert pool is None
                    yield None
                    return
                try:
                    (file_or_archive, destination_filename) = \
                        file_or_archive.split('#')
                except TypeError:
                    iface.fail(('%s is an invalid cache argument.'
                                    % file_or_archive),
                                steps=(job_flow[step_number:]
                                            if step_number != 0 else None))
                    failed = True
                    raise RuntimeError
                file_or_archive = os.path.expanduser(
                        os.path.expandvars(file_or_archive)
                    )
                file_or_archive_url = Url(file_or_archive)
                if not (file_or_archive_url.is_nfs
                            or file_or_archive_url.is_local):
                    iface.fail(('The file %s is not local or on NFS.'
                                    % file_or_archive),
                                steps=(job_flow[step_number:]
                                            if step_number != 0 else None))
                    failed = True
                    raise
                file_or_archive = file_or_archive_url.to_url()
                file_or_archive_basename = os.path.basename(file_or_archive)
                if archive:
                    archive_dir = destination_filename
                    destination_filename = os.path.basename(file_or_archive)
                if not os.path.isfile(file_or_archive):
                    iface.fail(('The file %s does not exist and thus cannot '
                                'be cached.') % file_or_archive,
                                steps=(job_flow[step_number:]
                                            if step_number != 0 else None))
                    failed = True
                    raise RuntimeError
                iface.status('Preparing temporary directories for storing '
                             '%s on slave nodes.' % file_or_archive_basename)
                '''Select engines to do "heavy lifting"; that is, they remove
                files copied to hosts on SIGINT/SIGTERM. Do it randomly
                (NO SEED) so if IWF occurs, second try will be different.
                IWF = intermittent weird failure. Set random seed so temp
                directory is reused if restarting job.'''
                random.seed(str(sorted(pid_map.keys())))
                engines_for_copying = [random.choice(list(engines)) 
                                        for engines in engine_map.values()
                                        if len(engines) > 0]
                '''Herd won't work with local engines; work around this by
                separating engines into two groups: local and remote.'''
                remote_hostnames_for_copying = list(
                        set(engine_map.keys()).difference(
                            set([current_hostname]))
                    )
                local_engines_for_copying = [
                        engine for engine in engines_for_copying
                        if engine in engine_map[current_hostname]
                    ]
                '''Create temporary directories on selected nodes; NOT
                WINDOWS-COMPATIBLE; must be changed if porting to Windows.'''
                if scratch == '-':
                    scratch_dir = tempfile.gettempdir()
                else:
                    scratch_dir = scratch
                temp_dir = os.path.join(
                                scratch_dir,
                                'dooplicity-%s' % ''.join(
                                        random.choice(string.ascii_uppercase
                                                        + string.digits)
                                        for _ in xrange(12)
                                    )
                            )
                '''To accommodate any slot-local BASH variables that may be in
                --scratch, echo them on all engines before adding to engine
                PYTHONPATHs.'''
                temp_dirs = apply_async_with_errors(pool, all_engines,
                    subprocess.check_output,
                    'echo "%s"' % temp_dir,
                    shell=True,
                    executable='/bin/bash',
                    message=('Error obtaining full paths of temporary '
                             'directories on cluster nodes. Restart IPython '
                             'engines and try again.'),
                    dict_format=True)
                for engine in temp_dirs:
                    temp_dirs[engine].strip()
                engines_with_unique_scratch, engines_to_symlink = [], []
                engine_to_copy_engine = {}
                for engine_for_copying in engines_for_copying:
                    for engine in engine_map[
                                host_map[engine_for_copying]
                            ]:
                        engine_to_copy_engine[engine] = engine_for_copying
                        if (engine != engine_for_copying
                            and temp_dirs[engine]
                            != temp_dirs[engine_for_copying]):
                            engines_with_unique_scratch.append(engine)
                            engines_to_symlink.append(engine)
                        elif engine == engine_for_copying:
                            engines_with_unique_scratch.append(engine)
                apply_async_with_errors(
                    pool, engines_for_copying, subprocess.check_output,
                    'mkdir -p %s' % temp_dir, shell=True,
                    executable='/bin/bash',
                    message=(('Error(s) encountered creating temporary '
                              'directories for storing {} on slave nodes. '
                              'Restart IPython engines and try again.').format(
                                                        file_or_archive
                                                    ))
                )
                if engines_to_symlink:
                    '''Create symlinks to resources in case of slot-local
                    scratch dirs'''
                    source_paths, destination_paths = {}, {}
                    for engine_to_symlink in engines_to_symlink:
                        source_paths[engine_to_symlink] = temp_dirs[
                                engine_to_copy_engine[engine_to_symlink]
                            ]
                        destination_paths[engine_to_symlink] = temp_dirs[
                                engine_to_symlink
                            ]
                    apply_async_with_errors(pool, engines_to_symlink,
                        os.remove, destination_paths,
                        message=('Error(s) encountered removing symlinks '
                                 'in slot-local scratch directories.'),
                        errors_to_ignore=['OSError'])
                    apply_async_with_errors(pool, engines_to_symlink,
                        os.symlink, source_paths, destination_paths,
                        message=('Error(s) encountered symlinking '
                                 'among slot-local scratch directories.'))
                # Add temp dirs to path
                apply_async_with_errors(
                    pool, all_engines, site.addsitedir, temp_dirs,
                    message=(('Error(s) encountered adding temporary '
                              'directories for storing {} to path on '
                              'slave nodes.').format(
                                                        file_or_archive
                                                    ))
                )
                '''Only foolproof way to die is by process polling. See
                http://stackoverflow.com/questions/284325/
                how-to-make-child-process-die-after-parent-exits for more
                information.'''
                apply_async_with_errors(
                    pool, engines_for_copying, subprocess.check_output,
                    ('echo "trap \\"{{ rm -rf {temp_dir}; exit 0; }}\\" '
                     'SIGHUP SIGINT SIGTERM EXIT; '
                     'while [[ \$(ps -p \$\$ -o ppid=) -gt 1 ]]; '
                     'do sleep 1; done & wait" '
                     '>{temp_dir}/delscript.sh').format(temp_dir=temp_dir),
                    shell=True,
                    executable='/bin/bash',
                    message=(
                            'Error creating script for scheduling temporary '
                            'directories on cluster nodes for deletion. '
                            'Restart IPython engines and try again.'
                        )
                )
                apply_async_with_errors(
                    pool, engines_for_copying, subprocess.Popen,
                    '/usr/bin/env bash %s/delscript.sh' % temp_dir, shell=True,
                    executable='/bin/bash',
                    message=(
                        'Error scheduling temporary directories on slave '
                        'nodes for deletion. Restart IPython engines and try '
                        'again.'
                    )
                )
                iface.status('Caching %s.' % file_or_archive_basename)
                destination_path = os.path.join(temp_dir, destination_filename)
                try:
                    import herd.herd as herd
                except ImportError:
                    '''Torrent distribution channel for compressed archive not
                    available.'''
                    apply_async_with_errors(
                        pool,
                        engines_for_copying,
                        shutil.copyfile,
                        file_or_archive, destination_path,
                        message=(('Error(s) encountered copying %s to '
                                  'slave nodes. Refer to the errors above '
                                  '-- and especially make sure $TMPDIR is not '
                                  'out of space on any node supporting an '
                                  'IPython engine -- before trying again.')
                                    % file_or_archive),
                    )
                else:
                    if local_engines_for_copying:
                        apply_async_with_errors(pool,
                            local_engines_for_copying,
                             subprocess.check_output, 'cp %s %s' % (
                                    file_or_archive, destination_path
                                ), shell=True, executable='/bin/bash',
                            message=(('Error(s) encountered copying %s to '
                                      'local filesystem. Refer to the errors '
                                      'above -- and especially make sure '
                                      '$TMPDIR is not out of space on any '
                                      'node supporting an IPython engine '
                                      '-- before trying again.')
                                        % file_or_archive),
                        )
                    if remote_hostnames_for_copying:
                        herd.run_with_opts(
                                file_or_archive,
                                destination_path,
                                hostlist=','.join(remote_hostnames_for_copying)
                            )
                # Extract if necessary
                if archive:
                    apply_async_with_errors(
                            pool, engines_for_copying, subprocess.check_output,
                            'rm -rf {}'.format(
                                    os.path.join(temp_dir, archive_dir)
                                ), shell=True, executable='/bin/bash'
                        )
                    apply_async_with_errors(
                            pool, engines_for_copying, subprocess.check_output,
                            'mkdir -p {}'.format(
                                os.path.join(temp_dir, archive_dir)
                            ), shell=True, executable='/bin/bash'
                    )
                    apply_async_with_errors(
                            pool, engines_for_copying, subprocess.check_output,
                            'tar xzf {} -C {}'.format(
                                destination_path,
                                os.path.join(temp_dir, archive_dir)),
                            shell=True,
                            executable='/bin/bash'
                    )
                    apply_async_with_errors(
                            pool, engines_for_copying, subprocess.check_output,
                            'rm -f {}'.format(destination_path),
                            shell=True, executable='/bin/bash'
                    )
                iface.step('Cached %s.' % file_or_archive_basename)
                try:
                    yield temp_dir
                finally:
                    # Cleanup
                    apply_async_with_errors(
                            pool,
                            engines_for_copying, subprocess.check_output,
                            'rm -rf {}'.format(temp_dir), shell=True,
                            executable='/bin/bash',
                            message=('Error(s) encountered removing temporary '
                                     'directories.')
                        )
        else:
            import multiprocessing
            def execute_balanced_job_with_retries(pool, iface,
                task_function, task_function_args,
                status_message='Tasks completed',
                finish_message='Completed tasks.', max_attempts=4):
                """ Executes parallel job locally with multiprocessing module.

                    Tasks are added to queue if they fail, and max_attempts-1
                    failures are permitted per task. 

                    pool: multiprocessing.Pool object
                    iface: DooplicityInterface object for spewing log messages
                        to console
                    task_function: name if function to execute
                    task_function_args: iterable of lists, each of whose
                        items are task_function's arguments, WITH THE EXCEPTION
                        OF A SINGLE KEYWORD ARGUMENT "attempt_count". This
                        argument must be the final keyword argument of the
                        function but _excluded_ from the arguments in any item
                        of task_function_args.
                    status_message: status message about tasks completed
                    finish_message: message to output when all tasks are
                        completed
                    max_attempts: max number of times to attempt any given
                        task

                    No return value.
                """
                global failed
                completed_tasks = 0
                tasks_to_assign = deque([
                        [task_function_arg, i, 0] for i, task_function_arg
                        in enumerate(task_function_args)
                    ])
                task_count = len(tasks_to_assign)
                assigned_tasks, asyncresults = {}, {}
                max_task_fails = 0
                iface.status(('    %s: %d/%d%s')
                                % (status_message, completed_tasks, task_count,
                                     (' | \\max_i (task_i fails): %d/%d'
                                       % (max_task_fails,
                                            max_attempts - 1)
                                       if max_attempts > 1 else '')))
                while completed_tasks < task_count:
                    if tasks_to_assign:
                        task_to_assign = tasks_to_assign.popleft()
                        asyncresults[task_to_assign[1]] = (
                                pool.apply_async(
                                    task_function,
                                    args=(task_to_assign[0] +
                                            [task_to_assign[2]])
                                )
                            )
                        assigned_tasks[task_to_assign[1]] = [
                                task_to_assign[0], task_to_assign[1],
                                task_to_assign[2] + 1
                            ]
                    asyncresults_to_remove = []
                    for task in asyncresults:
                        if asyncresults[task].ready():
                            return_value = asyncresults[task].get()
                            if return_value is not None:
                                if max_attempts > assigned_tasks[task][2]:
                                    # Add to queue for reattempt
                                    tasks_to_assign.append(
                                            assigned_tasks[task]
                                        )
                                    max_task_fails = max(
                                            assigned_tasks[task][2],
                                            max_task_fails
                                        )
                                    asyncresults_to_remove.append(task)
                                else:
                                    # Bail if max_attempts is saturated
                                    iface.fail(return_value,
                                    steps=(job_flow[step_number:]
                                            if step_number != 0 else None))
                                    failed = True
                                    raise RuntimeError
                            else:
                                # Success
                                completed_tasks += 1
                                asyncresults_to_remove.append(task)
                            iface.status(('    %s: %d/%d%s')
                                    % (status_message, completed_tasks,
                                        task_count,
                                        (' | \\max_i (task_i fails): %d/%d'
                                            % (max_task_fails,
                                                max_attempts - 1)
                                            if max_attempts > 1 else '')))
                    for task in asyncresults_to_remove:
                        del asyncresults[task]
                        del assigned_tasks[task]
                    time.sleep(0.1)
                iface.step(finish_message)
            @contextlib.contextmanager
            def cache(pool=None, file_or_archive=None, archive=True):
                """ Places X.[tar.gz/tgz]#Y in dir Y, unpacked if archive

                    pool: IPython Client object; all engines it spans are used
                    archive: file in format X.tar.gz#Y; False if nothing should
                        be done

                    Yields before deleting all files.
                """
                global failed
                if file_or_archive is None:
                    '''So with statements can be used all the time, even
                    when there's nothing to be archived'''
                    assert pool is None
                    yield None
                    return
                try:
                    (file_or_archive, destination_filename) = \
                        file_or_archive.split('#')
                except TypeError:
                    iface.fail(('%s is an invalid cache argument.'
                                    % file_or_archive),
                                steps=(job_flow[step_number:]
                                            if step_number != 0 else None))
                    failed = True
                    raise
                destination_filename = os.path.expanduser(
                        os.path.expandvars(destination_filename)
                    )
                file_or_archive_url = Url(file_or_archive)
                if not file_or_archive_url.is_local:
                    iface.fail(('The file %s is not local.'
                                    % file_or_archive),
                                steps=(job_flow[step_number:]
                                            if step_number != 0 else None))
                    failed = True
                    raise
                file_or_archive = file_or_archive_url.to_url()
                file_or_archive_basename = os.path.basename(file_or_archive)
                if archive:
                    archive_dir = destination_filename
                    destination_filename = os.path.basename(file_or_archive)
                if not os.path.isfile(file_or_archive):
                    iface.fail(('The file %s does not exist and thus cannot '
                                'be cached.') % file_or_archive,
                                steps=(job_flow[step_number:]
                                            if step_number != 0 else None))
                    failed = True
                    raise RuntimeError
                temp_dir = make_temp_dir_and_register_cleanup(
                                                    None if scratch == '-'
                                                    else scratch
                                                )
                iface.status('Caching %s.' % file_or_archive_basename)
                destination_path = os.path.join(temp_dir, destination_filename)
                shutil.copyfile(file_or_archive, destination_path)
                # Extract if necessary
                if archive:
                    try:
                        os.makedirs(os.path.join(temp_dir, archive_dir))
                    except OSError:
                        # Hopefully, directory is already created
                        pass
                    import subprocess
                    try:
                        subprocess.check_output(
                                        'tar xzf {} -C {}'.format(
                                            destination_path,
                                            os.path.join(temp_dir, archive_dir)
                                        ),
                                        shell=True,
                                        bufsize=-1,
                                        stderr=subprocess.STDOUT,
                                        executable='/bin/bash'
                                    )
                    except subprocess.CalledProcessError as e:
                        iface.fail(
                                ('Decompression of archive failed; exit code '
                                 'was %s, and reason was "%s".') % (
                                    e.returncode,
                                    e.output.strip()
                                ),
                                steps=(job_flow[step_number:]
                                            if step_number != 0 else None)
                            )
                        failed = True
                        raise RuntimeError
                    try:
                        os.remove(destination_path)
                    except OSError:
                        pass
                iface.step('Cached %s.' % file_or_archive_basename)
                try:
                    yield temp_dir
                finally:
                    # Cleanup
                    shutil.rmtree(temp_dir)
        # Serialize JSON configuration
        if json_config is not None:
            with open(json_config) as json_stream:
                full_payload = json.load(json_stream)
        else:
            full_payload = json.load(sys.stdin)
        try:
            job_flow = full_payload['Steps']
        except KeyError:
            iface.fail(
                    'Input JSON not in proper format. Ensure that the JSON '
                    'object has a Steps key.'
                )
            failed = True
            raise
        step_count = len(job_flow)
        steps = OrderedDict()
        try:
            for step in job_flow:
                step_args = {}
                j = 0
                j_max = len(step['HadoopJarStep']['Args'])
                while j < j_max:
                    arg_name = step['HadoopJarStep']['Args'][j][1:].strip()
                    if arg_name == 'D':
                        D_arg = step['HadoopJarStep']['Args'][j+1].split('=')
                        if D_arg[0] in ['mapred.reduce.tasks',
                                        'mapreduce.job.reduces']:
                            step_args['task_count'] = int(D_arg[1])
                        elif D_arg[0] \
                            in ['mapred.text.key.partitioner.options',
                                'mapreduce.partition.keypartitioner.options']:
                            step_args['partition_options'] = D_arg[1]
                        elif D_arg[0] \
                            == 'stream.num.map.output.key.fields':
                            step_args['key_fields'] = int(D_arg[1])
                            # Default sort is now across key fields
                            if 'sort_options' not in step_args:
                                step_args['sort_options'] = (
                                            '-k1,%d' % step_args['key_fields']
                                        )
                        elif D_arg[0] \
                            in ['mapred.text.key.comparator.options',
                                'mapreduce.partition.keycomparator.options']:
                            step_args['sort_options'] = D_arg[1]
                        j += 2
                    elif arg_name == 'input':
                        try:
                            step_args['input'] = ','.join(
                                    [step['HadoopJarStep']['Args'][j+1],
                                     step_args['input']]
                                )
                        except KeyError:
                            step_args['input'] \
                                = step['HadoopJarStep']['Args'][j+1]
                        j += 2
                    elif arg_name == 'multiOutput':
                        step_args['multiple_outputs'] = True
                        j += 1
                    elif arg_name == 'lazyOutput':
                        # Do nothing
                        j += 1
                    else:
                        step_args[step['HadoopJarStep']['Args'][j][1:]] \
                            = step['HadoopJarStep']['Args'][j+1].strip()
                        j += 2
                # Set default options
                if 'key_fields' not in step_args:
                    step_args['key_fields'] = 1
                if 'partition_options' not in step_args:
                    step_args['partition_options'] = '-k1'
                if 'sort_options' not in step_args:
                    step_args['sort_options'] = '-k1'
                steps[step['Name']] = step_args
        except (KeyError, IndexError):
            iface.fail(
                    'JSON file not in proper format. Ensure '
                    'that each step object has a HadoopJarStep '
                    'object with an Args array and a Name string.'
                )
            failed = True
            raise
        '''Check steps for required Hadoop streaming command-line parameters
        and for whether outputs are writable.'''
        missing_data = defaultdict(list)
        bad_output_data = []
        required_data = set(['input', 'output', 'mapper', 'reducer'])
        identity_mappers \
            = set(['cat', 'org.apache.hadoop.mapred.lib.IdentityMapper'])
        identity_reducers \
            = set(['cat', 'org.apache.hadoop.mapred.lib.IdentityReducer'])
        errors = []
        for step in steps:
            step_data = steps[step]
            for required_parameter in required_data:
                if required_parameter not in step_data:
                    missing_data[step].append('-' + required_parameter)
                elif not force and required_parameter == 'output' \
                    and os.path.exists(step_data['output']):
                    bad_output_data.append(step)
            try:
                if step_data['inputformat'] \
                    == 'org.apache.hadoop.mapred.lib.NLineInputFormat' \
                    and os.path.isdir(step_data['input']):
                    errors.append(('In step "%s", input should be a single '
                                   'file if using NLineFormat, but '
                                   '"%s" was specified.') % (
                                                        step,
                                                        step_data['input']
                                                    ))
            except KeyError:
                pass
        if missing_data:
            errors.extend(['Step "%s" is missing required parameter(s) "%s".' % 
                                (step, ', '.join(missing_data[step]))
                                for step in missing_data])
        if bad_output_data:
            errors.extend(['Output directory name "%s" of step "%s" already '
                           'exists as a file or directory, and --force was '
                           'not invoked to permit overwriting it.' 
                           % (steps[step]['output'], step)
                           for step in bad_output_data])
        if errors:
            iface.fail('\n'.join([(('%d) ' % (i+1)) + error)
                                    if len(errors) > 1 else errors[0]
                                    for i, error in enumerate(errors)]))
            failed = True
            raise RuntimeError
        if not keep_intermediates:
            # Create schedule for deleting intermediates
            marked_intermediates = set()
            all_outputs = set()
            post_step_cleanups = defaultdict(list)
            '''Traverse steps in reverse order to obtain when an intermediate
            directory is last used.'''
            for i, step in enumerate(
                                OrderedDict(reversed(steps.items()[1:]))
                            ):
                step_inputs = [os.path.abspath(step_input) for step_input in
                                steps[step]['input'].split(',')]
                all_outputs.add(os.path.abspath(steps[step]['output']))
                for step_input in step_inputs:
                    if step_input not in marked_intermediates:
                        post_step_cleanups[step_count - i - 1].append(
                                step_input
                            )
                        marked_intermediates.add(step_input)
        # Create intermediate directories
        for step in steps:
            try:
                shutil.rmtree(steps[step]['output'])
            except OSError:
                # May be a file then
                try:
                    os.remove(steps[step]['output'])
                except OSError:
                    # Just didn't exist
                    pass
            try:
                os.makedirs(steps[step]['output'])
            except OSError:
                iface.fail(('Problem encountered trying to create '
                            'directory %s.') % steps[step]['output'])
                failed = True
                raise
            map_err_dir = os.path.join(steps[step]['output'], 'dp.map.log')
            try:
                os.makedirs(map_err_dir)
            except OSError:
                iface.fail(('Problem encountered trying to create '
                            'directory %s.') % map_err_dir)
                failed = True
                raise
            reduce_err_dir = os.path.join(steps[step]['output'],
                                            'dp.reduce.log')
            try:
                os.makedirs(reduce_err_dir)
            except OSError:
                iface.fail(('Problem encountered trying to create '
                            'directory %s.') % reduce_err_dir)
                failed = True
                raise
        # Run steps
        step_number = 0
        total_steps = len(steps)
        if not ipy:
            # Pool's only for if we're in local mode
            try:
                pool = multiprocessing.Pool(num_processes, init_worker,
                                                maxtasksperchild=5)
            except Exception:
                # maxtasksperchild doesn't work, somehow? Supported only in 2.7
                pool = multiprocessing.Pool(num_processes, init_worker)
        for step in steps:
            step_data = steps[step]
            step_inputs = []
            # Handle multiple input files/directories
            for input_file_or_dir in step_data['input'].split(','):
                if os.path.isfile(input_file_or_dir):
                    step_inputs.append(input_file_or_dir)
                elif os.path.isdir(input_file_or_dir):
                    step_inputs.extend(
                            glob.glob(os.path.join(input_file_or_dir, '*'))
                        )
            # TODO: support cacheArchives and cacheFile simultaneously
            if 'archives' in step_data or 'cacheArchive' in step_data:
                # Prefer archives to cacheArchives
                try:
                    to_cache = step_data['archives']
                except KeyError:
                    to_cache = step_data['cacheArchive']
            elif 'files' in step_data or 'cacheFile' in step_data:
                try:
                    to_cache = step_data['files']
                except KeyError:
                    to_cache = step_data['cacheFile']
            else:
                to_cache = None
            with cache(pool if to_cache else None, to_cache,
                        True if 'archives'
                        in step_data else False) as dir_to_path:
                if step_data['mapper'] not in identity_mappers:
                    # Perform map step only if mapper isn't identity
                    try:
                        if ('multiple_outputs' in step_data) or \
                            step_data['outputformat'] \
                            in ['edu.jhu.cs.MultipleOutputFormat',
                                'edu.jhu.cs.'
                                'MultipleIndexedLzoTextOutputFormat']:
                            multiple_outputs = True
                        else:
                            multiple_outputs = False
                    except KeyError:
                        # No multiple outputs
                        multiple_outputs = False
                    try:
                        if step_data['inputformat'] \
                            == 'org.apache.hadoop.mapred.lib.NLineInputFormat':
                            nline_input = True
                        else:
                            nline_input = False
                    except KeyError:
                        # Don't assign one line per mapper
                        nline_input = False
                    output_dir = step_data['output']
                    if step_data['reducer'] not in identity_reducers:
                        '''There's a reducer parameter, so input to reducer is
                        output of mapper. Change output directory.'''
                        output_dir = os.path.join(output_dir, 'dp.map')
                        try:
                            os.makedirs(output_dir)
                        except OSError:
                            if os.path.exists(output_dir):
                                pass
                            else:
                                iface.fail(('Problem encountered trying to '
                                            'create directory %s.')
                                            % output_dir,
                                            steps=(job_flow[step_number:]
                                                       if step_number != 0
                                                       else None))
                                failed = True
                                raise
                        try:
                            if ('multiple_outputs' in step_data) or \
                                step_data['outputformat'] \
                                in ['edu.jhu.cs.MultipleOutputFormat',
                                    'edu.jhu.cs.'
                                    'MultipleIndexedLzoTextOutputFormat']:
                                # Multiple outputs apply AFTER reduce step
                                multiple_outputs = False
                        except KeyError:
                            # No outputformat
                            pass
                    if nline_input:
                        # Create temporary input files
                        split_input_dir = make_temp_dir(common)
                        input_files = []
                        try:
                            with open(step_inputs[0]) as nline_stream:
                                for i, line in enumerate(nline_stream):
                                    offset = str(i)
                                    input_files.append(os.path.join(
                                                            split_input_dir,
                                                            offset
                                                        )
                                                    )
                                    with open(input_files[-1], 'w') \
                                        as output_stream:
                                        print >>output_stream, separator.join([
                                                                offset, line
                                                            ])
                        except IndexError:
                            raise RuntimeError('No NLineInputFormat input to '
                                               'step "%s".' % step)
                    else:
                        input_files = [input_file for input_file in step_inputs
                                        if os.path.isfile(input_file)]
                    input_file_count = len(input_files)
                    if not input_file_count:
                        iface.step('No input found; skipping step.')
                    err_dir = os.path.join(steps[step]['output'], 'dp.map.log')
                    iface.step('Step %d/%d: %s' % 
                                (step_number + 1, total_steps, step))
                    iface.status('    Starting step runner...')
                    execute_balanced_job_with_retries(
                            pool, iface, step_runner_with_error_return,
                                       [[step_data['mapper'], input_file,
                                         output_dir, err_dir,
                                         i, multiple_outputs,
                                         separator, None, None, gzip,
                                         gzip_level, scratch, direct_write,
                                         sort, dir_to_path]
                                         for i, input_file
                                         in enumerate(input_files)
                                         if os.path.isfile(input_file)],
                            status_message='Tasks completed',
                            finish_message=(
                                '    Completed %s.'
                                % dp_iface.inflected(input_file_count, 'task')
                            ),
                            max_attempts=max_attempts
                        )
                    # Adjust step inputs in case a reducer follows
                    step_inputs = [input_file for input_file 
                                    in glob.glob(output_dir)
                                    if os.path.isfile(input_file)]
                if step_data['reducer'] not in identity_reducers:
                    '''Determine whether to use "mod" partitioner that uses
                    product of key fields % reducer count to assign tasks.'''
                    try:
                        if (step_data['partitioner']
                                == 'edu.jhu.cs.ModPartitioner'):
                            mod_partition = True
                        else:
                            mod_partition = False
                    except KeyError:
                        # Default to no mod partition
                        mod_partition = False
                    # Partition inputs into tasks, presorting
                    output_dir = os.path.join(step_data['output'], 'dp.tasks')
                    try:
                        os.makedirs(output_dir)
                    except OSError:
                        if os.path.exists(output_dir):
                            pass
                        else:
                            iface.fail(('Problem encountered trying to '
                                        'create directory %s.') % output_dir,
                                        steps=(job_flow[step_number:]
                                                if step_number != 0 else None))
                            failed = True
                            raise
                    input_files = [input_file for input_file in step_inputs
                                    if os.path.isfile(input_file)]
                    input_file_count = len(input_files)
                    if input_file_count > num_processes:
                        file_count_per_group = input_file_count / num_processes
                        input_file_groups = [
                                input_files[k:k+file_count_per_group]
                                for k in
                                xrange(0, input_file_count,
                                        file_count_per_group)
                            ]
                    else:
                        input_file_groups = [[input_file]
                                                for input_file in input_files]
                    input_file_group_count = len(input_file_groups)
                    iface.step('Step %d/%d: %s'
                                 % (step_number + 1, total_steps, step))
                    execute_balanced_job_with_retries(
                            pool, iface, presorted_tasks,
                            [[input_file_group, i,
                                step_data['sort_options'], output_dir,
                                step_data['key_fields'], separator,
                                step_data['partition_options'],
                                step_data['task_count'], memcap, gzip,
                                gzip_level, scratch, direct_write,
                                sort, mod_partition]
                                    for i, input_file_group
                                    in enumerate(input_file_groups)],
                            status_message='Inputs partitioned',
                            finish_message=(
                                '    Partitioned %s into tasks.'
                                % dp_iface.inflected(input_file_group_count,
                                                     'input')
                            ),
                            max_attempts=max_attempts
                        )
                    iface.status('    Starting step runner...')
                    input_files = [os.path.join(output_dir, '%d.*' % i) 
                                   for i in xrange(step_data['task_count'])]
                    # Filter out bad globs
                    input_files = [input_file for input_file in input_files
                                    if glob.glob(input_file)]
                    input_file_count = len(input_files)
                    try:
                        multiple_outputs = (
                                ('multiple_outputs' in step_data) or
                                    step_data['outputformat']
                                    in ['edu.jhu.cs.MultipleOutputFormat',
                                        'edu.jhu.cs.'
                                        'MultipleIndexedLzoTextOutputFormat']
                            )
                    except KeyError:
                        multiple_outputs = False
                    err_dir = os.path.join(
                                    steps[step]['output'],
                                    'dp.reduce.log'
                                )
                    output_dir = step_data['output']
                    return_values = []
                    execute_balanced_job_with_retries(
                            pool, iface, step_runner_with_error_return,
                                [[step_data['reducer'], input_file, output_dir, 
                                err_dir, i, multiple_outputs, separator,
                                step_data['sort_options'], memcap, gzip,
                                gzip_level, scratch, direct_write,
                                sort, dir_to_path]
                                    for i, input_file
                                    in enumerate(input_files)],
                            status_message='Tasks completed',
                            finish_message=(
                                '    Completed %s.'
                                % dp_iface.inflected(input_file_count, 'task')
                            ),
                            max_attempts=max_attempts
                        )
            # Really close open file handles in PyPy
            gc.collect()
            if not keep_intermediates:
                iface.status('    Deleting temporary files...')
                # Kill NLineInput files if they're there
                try:
                    shutil.rmtree(split_input_dir)
                except (NameError, OSError):
                    pass
                try:
                    # Intermediate map output should be deleted if it exists
                    shutil.rmtree(
                                os.path.join(
                                        step_data['output'], 'dp.map'
                                    )
                            )
                except OSError:
                    pass
                try:
                    # Remove dp.tasks directory
                    shutil.rmtree(
                            os.path.join(step_data['output'], 'dp.tasks')
                        )
                except OSError:
                    pass
                for to_remove in post_step_cleanups[step_number]:
                    if to_remove not in all_outputs:
                        '''Remove directory only if it's an -output of some
                        step and an -input of another step.'''
                        continue
                    if os.path.isfile(to_remove):
                        try:
                            os.remove(to_remove)
                        except OSError:
                            pass
                    elif os.path.isdir(to_remove):
                        for detritus in glob.iglob(
                                            os.path.join(to_remove, '*')
                                        ):
                            if detritus[-4:] != '.log':
                                try:
                                    os.remove(detritus)
                                except OSError:
                                    try:
                                        shutil.rmtree(detritus)
                                    except OSError:
                                        pass
                        if not os.listdir(to_remove):
                            try:
                                os.rmdir(to_remove)
                            except OSError:
                                pass
                iface.step('    Deleted temporary files.')
            step_number += 1
        if not ipy:
            pool.close()
        if not keep_last_output and not keep_intermediates:
            try:
                os.remove(step_data['output'])
            except OSError:
                # Not a file; treat as dir
                for detritus in glob.iglob(
                                os.path.join(step_data['output'], '*')
                            ):
                    if detritus[-4:] != '.log':
                        try:
                            os.remove(detritus)
                        except OSError:
                            # Not a file
                            try:
                                shutil.rmtree(detritus)
                            except OSError:
                                # Phantom; maybe user deleted it
                                pass
        if not keep_intermediates:
            for step in steps:
                step_data = steps[step]
                try:
                    os.remove(step_data['output'])
                except OSError:
                    # Not a file; treat as dir
                    for detritus in glob.iglob(
                                    os.path.join(step_data['output'], '*')
                                ):
                        if detritus[-4:] != '.log':
                            try:
                                os.remove(detritus)
                            except OSError:
                                # Not a file
                                try:
                                    shutil.rmtree(detritus)
                                except OSError:
                                    # Phantom; maybe user deleted it
                                    pass
        iface.done()
    except (Exception, GeneratorExit):
        # GeneratorExit added just in case this happens on modifying code
        if 'interrupt_engines' in locals():
            interrupt_engines(pool, iface)
        if not failed:
            time.sleep(0.2)
            if 'step_number' in locals():
                iface.fail(steps=(job_flow[step_number:]
                            if step_number != 0 else None))
            else:
                iface.fail()
        if 'split_input_dir' in locals():
            '''raise below refers to last exception, so can't try-except
            OSError here'''
            if os.path.isdir(split_input_dir):
                shutil.rmtree(split_input_dir)
        raise
    except (KeyboardInterrupt, SystemExit):
        if 'interrupt_engines' in locals():
            interrupt_engines(pool, iface)
        if 'step_number' in locals():
            iface.fail(steps=(job_flow[step_number:]
                        if step_number != 0 else None),
                        opener='*****Terminated*****')
        else:
            iface.fail()
        if 'pool' in locals() and 'interrupt_engines' not in locals():
            pool.terminate()
            pool.join()
        if 'split_input_dir' in locals():
            try:
                shutil.rmtree(split_input_dir)
            except OSError:
                pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, 
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    add_args(parser)
    dp_iface.add_args(parser)
    args = parser.parse_args(sys.argv[1:])
    run_simulation(args.branding, args.json_config, args.force,
                    args.memcap, args.num_processes, args.separator,
                    args.keep_intermediates, args.keep_last_output,
                    args.log, args.gzip_outputs, args.gzip_level,
                    args.ipy, args.ipcontroller_json, args.ipy_profile,
                    args.scratch, args.common, args.sort, args.max_attempts,
                    args.direct_write)