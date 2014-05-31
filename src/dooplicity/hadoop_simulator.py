#!/usr/bin/env python
"""
hadoop_simulator.py
Part of Dooplicity framework

Runs JSON-encoded Hadoop Streaming job flow. FUNCTIONALITY IS IDIOSYNCRATIC;
it is currently confined to those features used by Rail. Format of input JSON
mirrors that of StepConfig list from JSON sent to EMR via RunJobsFlow. Any
files input to a mapper can be gzip'd, but inputs to a reducer currently cannot
be.

TODO: Implement -k option: Kill intermediates when they're no longer necessary.
This requires checking rest of steps to see if an intermediate directory
appears.

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
import os
import shutil
from collections import defaultdict, OrderedDict
import time
import glob
import hashlib
import subprocess
import json
import interface as dp_iface
import signal
import gc

def init_worker():
    """ Prevents KeyboardInterrupt from reaching a pool's workers.

        Exiting gracefully after KeyboardInterrupt or SystemExit is a
        challenge. The solution implemented here is by John Reese and is from
        http://noswap.com/blog/python-multiprocessing-keyboardinterrupt .

        No return value.
    """
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def presorted_tasks(input_files, process_id, sort_options, output_dir,
                    key_fields, separator, task_count, memcap):
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
        streaming_command: streaming command to run.
        task_count: number of tasks in which to partition input.
        memcap: maximum amount of memory to use per UNIX sort instance;
            ignored if sort is None.

        Return value: Empty tuple if no errors encountered; otherwise
            (input_files,).
    """
    try:
        task_streams = {}
        for input_file in input_files:
            with open(input_file) as input_stream:
                for line in input_stream:
                    key = separator.join(
                                line.strip().split(separator)[:key_fields]
                            )
                    task = int(hashlib.md5(key).hexdigest(), 16) % task_count
                    try:
                        task_streams[task].write(line)
                    except KeyError:
                        # Task file doesn't exist yet; create it
                        task_file = os.path.join(output_dir, str(task) + '.'
                                                             + str(process_id)
                                                             + '.unsorted')
                        task_streams[task] = open(task_file, 'w')
                        task_streams[task].write(line)
        for task in task_streams:
            task_streams[task].close()
        # Presort task files
        for unsorted_file in glob.glob(os.path.join(
                                                output_dir,
                                                '*.%s.unsorted' % process_id
                                            )):
            sort_command = 'sort -S %d %s %s > %s' % (memcap,
                                                        sort_options,
                                                        unsorted_file,
                                                        unsorted_file[:-9])
            sort_return = subprocess.call(sort_command, shell=True, bufsize=-1)
            if sort_return != 0:
                return ('Error encountered sorting input file %s.' %
                            input_file,)
            os.remove(unsorted_file)
        return tuple()
    except Exception as e:
        # Uncaught miscellaneous exception
        return ('Error "%s" encountered partitioning input files '
                '[%s] into tasks.'
                        % (e.message, (('%s, '* (len(input_files) - 1) 
                                       + '%s') % tuple(input_files))),)

def step_runner_with_error_return(streaming_command, input_glob, output_dir,
                                  err_dir, task_id, multiple_outputs,
                                  separator, sort_options, memcap):
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
        memcap: maximum amount of memory to use per UNIX sort instance;
            ignored if sort is None.

        Return value: empty tuple iff step runs successfully; otherwise, tuple
            (error message, full streaming command that was run)
    """
    command_to_run = None
    try:
        input_files = [input_file for input_file in glob.glob(input_glob)
                        if os.path.isfile(input_file)]
        if not input_files:
            # No input!
            return tuple()
        if sort_options is None:
            # Mapper. Check if first input file is gzip'd
            with open(input_files[0], 'rb') as binary_input_stream:
                if binary_input_stream.read(2) == '\x1f\x8b':
                    # Magic number of gzip'd file found
                    prefix = 'gunzip -c %s' % input_glob
                else:
                    prefix = 'cat %s' % input_glob
        else:
            # Reducer. Merge sort the input glob.
            prefix = 'sort -S %d %s -m %s' % (memcap, sort_options, input_glob)
        err_file = os.path.join(err_dir, '%d.log' % task_id)
        new_env = os.environ.copy()
        new_env['mapred_task_partition'] = str(task_id)
        if multiple_outputs:
            # Must grab each line of output and separate by directory
            command_to_run \
                = prefix + ' | ' + streaming_command + (' 2>%s' % err_file)
            multiple_output_process = subprocess.Popen(
                    command_to_run,
                    shell=True,
                    stdout=subprocess.PIPE,
                    env=new_env,
                    bufsize=-1
                )
            task_file_streams = {}
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
                        if os.path.exists(key_dir):
                            pass
                        else:
                            return (('Problem encountered trying to create '
                                     'directory %s.') % key_dir,
                                    command_to_run)
                    task_file_streams[key] = open(os.path.join(key_dir,
                                                        str(task_id)), 'w')
                    task_file_streams[key].write(line_to_write)
            multiple_output_process_return = multiple_output_process.wait()
            if multiple_output_process_return != 0:
                return (('Exit level was %d.')
                        % multiple_output_process_return, command_to_run)
            for key in task_file_streams:
                task_file_streams[key].close()
        else:
            out_file = os.path.join(output_dir, str(task_id))
            command_to_run \
                = prefix + ' | ' + streaming_command + (' >%s 2>%s'
                                                             % (out_file,
                                                                err_file))
            single_output_process_return = subprocess.call(
                                                    command_to_run,
                                                    shell=True,
                                                    env=new_env,
                                                    bufsize=-1
                                                )
            if single_output_process_return != 0:
                return (('Exit level was %d.')
                        % single_output_process_return, command_to_run)
        return tuple()
    except Exception as e:
        # Uncaught miscellaneous exception
        return ('Error "%s" executing task on input %s'
                % (e.message, input_file), command_to_run)

def run_simulation(branding, json_config, force, memcap, num_processes,
                    separator, keep_intermediates, keep_last_output,
                    log):
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
        memcap: maximum amount of memory to use per UNIX sort instance.
        num_processes: number of subprocesses to open at once.
        separator: separator between successive fields in inputs and 
            intermediates.
        keep_intermediates: keeps all intermediate output.
        keep_last_output: keeps outputs that are unused as inputs by steps.
        log: name of file in which to store messages written to stderr

        No return value.
    """
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
                for j in xrange(0, len(step['HadoopJarStep']['Args']), 2):
                    arg_name = step['HadoopJarStep']['Args'][j][1:].strip()
                    if arg_name == 'D':
                        D_arg = step['HadoopJarStep']['Args'][j+1].split('=')
                        if D_arg[0] == 'mapred.reduce.tasks':
                            step_args['task_count'] = int(D_arg[1])
                        elif D_arg[0] == 'mapred.text.key.partitioner.options':
                            step_args['key_fields'] \
                                = int(D_arg[1].split(',')[-1])
                        elif D_arg[0] == 'stream.num.map.output.key.fields':
                            step_args['sort_options'] = '-k1,%s' % D_arg[1]
                    elif arg_name == 'input':
                        try:
                            step_args['input'] = ','.join(
                                    [step['HadoopJarStep']['Args'][j+1],
                                     step_args['input']]
                                )
                        except KeyError:
                            step_args['input'] \
                                = step['HadoopJarStep']['Args'][j+1]
                    else:
                        step_args[step['HadoopJarStep']['Args'][j][1:]] \
                            = step['HadoopJarStep']['Args'][j+1].strip()
                steps[step['Name']] = step_args
        except KeyError, IndexError:
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
        for step in steps:
            step_data = steps[step]
            for required_parameter in required_data:
                if required_parameter not in step_data:
                    missing_data[step].append('-' + required_parameter)
                elif not force and required_parameter == 'output' \
                    and os.path.exists(step_data['output']):
                    bad_output_data.append(step)
        errors = []
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
            iface.fail('\n'.join([('%d) ' % (i+1)) + error
                                    for i, error in enumerate(errors)]))
            failed = True
            raise RuntimeError
        if not keep_intermediates:
            # Create schedule for deleting intermediates
            marked_intermediates = set()
            post_step_cleanups = defaultdict(list)
            '''Traverse steps in reverse order to obtain when an intermediate
            directory is last used. DON'T MESS WITH FIRST INPUT.'''
            for i, step in enumerate(
                                OrderedDict(reversed(steps.items()[1:]))
                            ):
                step_inputs = [os.path.abspath(step_input) for step_input in
                                steps[step]['input'].split(',')]
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
        import multiprocessing
        step_number = 0
        total_steps = len(steps)
        for step in steps:
            step_data = steps[step]
            step_inputs = []
            # Handle multiple input files/directories
            for input_file_or_dir in step_data['input'].split(','):
                if os.path.isfile(input_file_or_dir):
                    step_inputs.append(input_file_or_dir)
                else:
                    step_inputs.extend(
                            glob.glob(os.path.join(input_file_or_dir, '*'))
                        )
            if step_data['mapper'] not in identity_mappers:
                # Perform map step only if mapper isn't identity
                try:
                    if step_data['outputformat'] \
                        == 'edu.jhu.cs.MultipleOutputFormat':
                        multiple_outputs = True
                    else:
                        multiple_outputs = False
                except KeyError:
                    # No multiple outputs
                    multiple_outputs = False
                    pass
                output_dir = step_data['output']
                if step_data['reducer'] not in identity_mappers:
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
                                        'create directory %s.') % output_dir,
                                        steps=(job_flow[step_number:]
                                                   if step_number != 0
                                                   else None))
                            failed = True
                            raise
                    try:
                        if step_data['outputformat'] \
                            == 'edu.jhu.cs.MultipleOutputFormat':
                            # Multiple outputs apply AFTER reduce step
                            multiple_outputs = False
                    except KeyError:
                        # No outputformat
                        pass
                return_values = []
                pool = multiprocessing.Pool(num_processes, init_worker)
                input_files = [input_file for input_file in step_inputs
                                if os.path.isfile(input_file)]
                input_file_count = len(input_files)
                err_dir = os.path.join(steps[step]['output'], 'dp.map.log')
                for i, input_file in enumerate(input_files):
                    if os.path.isfile(input_file):
                        pool.apply_async(
                            step_runner_with_error_return, 
                            args=(step_data['mapper'], input_file,
                                output_dir, err_dir, i,
                                multiple_outputs, 
                                separator,
                                None,
                                None),
                            callback=return_values.append
                        )
                pool.close()
                iface.step('Step %d/%d: %s' % 
                            (step_number + 1, total_steps, step))
                while len(return_values) != input_file_count:
                    try:
                        max_tuple = max(map(len, return_values))
                    except ValueError:
                        # return_values is empty
                        max_tuple = -1
                        pass
                    if max_tuple > 0:
                        '''There are error tuples; could put error message
                        directly in RuntimeError, but then it would be
                        positioned after the Dooplicity message about
                        resuming the job.'''
                        errors = ['Streaming command "%s" failed: %s' % 
                                  (error[1], error[0]) for error
                                  in return_values if len(error) == 2]
                        errors = [('%d) ' % (i + 1)) + error
                                    for i, error in enumerate(errors)]
                        iface.fail('\n'.join(errors),
                                    steps=(job_flow[step_number:]
                                            if step_number != 0 else None))
                        failed = True
                        raise RuntimeError
                    iface.status('    Tasks completed: %d/%d'
                                 % (len(return_values), input_file_count))
                    time.sleep(0.2)
                try:
                    max_tuple = max(map(len, return_values))
                except ValueError:
                    # return_values is empty
                    max_tuple = -1
                    pass
                if max_tuple > 0:
                    '''There are error tuples; could put error message
                    directly in RuntimeError, but then it would be
                    positioned after the Dooplicity message about
                    resuming the job.'''
                    errors = ['Streaming command "%s" failed: %s' % 
                              (error[1], error[0]) for error
                              in return_values if len(error) == 2]
                    errors = [('%d) ' % (i + 1)) + error
                                for i, error in enumerate(errors)]
                    iface.fail('\n'.join(errors),
                                steps=(job_flow[step_number:]
                                        if step_number != 0 else None))
                    failed = True
                    raise RuntimeError
                iface.step('    Completed %s.'
                           % dp_iface.inflected(input_file_count, 'task'))
                # Adjust step inputs in case a reducer follows
                step_inputs = [input_file for input_file 
                                in glob.glob(output_dir)
                                if os.path.isfile(input_file)]
            if step_data['reducer'] not in identity_mappers:
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
                return_values = []
                input_files = [input_file for input_file in step_inputs
                                if os.path.isfile(input_file)]
                input_file_count = len(input_files)
                if input_file_count > num_processes:
                    file_count_per_group = input_file_count / num_processes
                    input_file_groups = [
                            input_files[k:k+file_count_per_group]
                            for k in
                            xrange(0, input_file_count, file_count_per_group)
                        ]
                else:
                    input_file_groups = [[input_file]
                                            for input_file in input_files]
                input_file_group_count = len(input_file_groups)
                pool = multiprocessing.Pool(num_processes, init_worker)
                for i, input_file_group in enumerate(input_file_groups):
                    pool.apply_async(
                                        presorted_tasks,
                                        args=(input_file_group, i,
                                            step_data['sort_options'],
                                            output_dir,
                                            step_data['key_fields'],
                                            separator,
                                            step_data['task_count'],
                                            memcap),
                                        callback=return_values.append
                                    )
                pool.close()
                iface.step('Step %d/%d: %s'
                             % (step_number + 1, total_steps, step))
                while len(return_values) != input_file_group_count:
                    try:
                        max_tuple = max(map(len, return_values))
                    except ValueError:
                        # return_values is empty
                        max_tuple = -1
                        pass
                    if max_tuple > 0:
                        # There are error tuples
                        errors = [error[0] for error in return_values if error]
                        errors = [('%d) ' % (i + 1)) + error
                                    for i, error in enumerate(errors)]
                        iface.fail('\n'.join(errors),
                                   (job_flow[step_number:]
                                    if step_number != 0 else None))
                        failed = True
                        raise RuntimeError
                    iface.status(('    Inputs partitioned: %d/%d\r')
                                 % (len(return_values),
                                        input_file_group_count))
                    time.sleep(0.2)
                try:
                    max_tuple = max(map(len, return_values))
                except ValueError:
                    # return_values is empty
                    max_tuple = -1
                    pass
                if max_tuple > 0:
                    # There are error tuples
                    errors = [error[0] for error in return_values if error]
                    errors = [('%d) ' % (i + 1)) + error
                                for i, error in enumerate(errors)]
                    iface.fail('\n'.join(errors),
                               (job_flow[step_number:]
                                if step_number != 0 else None))
                    failed = True
                    raise RuntimeError
                iface.step('    Partitioned %s into tasks.'
                            % dp_iface.inflected(input_file_group_count,
                                                 'input'))
                return_values = []
                input_files = [os.path.join(output_dir, '%d.*' % i) 
                               for i in xrange(step_data['task_count'])]
                # Filter out bad globs
                input_files = [input_file for input_file in input_files
                                if glob.glob(input_file)]
                input_file_count = len(input_files)
                pool = multiprocessing.Pool(num_processes, init_worker)
                try:
                    multiple_outputs = (step_data['outputformat']
                                        == 'edu.jhu.cs.MultipleOutputFormat')
                except KeyError:
                    multiple_outputs = False
                err_dir = os.path.join(steps[step]['output'], 'dp.reduce.log')
                output_dir = step_data['output']
                for i, input_file in enumerate(input_files):
                    pool.apply_async(
                            step_runner_with_error_return, 
                            args=(step_data['reducer'], input_file,
                                output_dir, err_dir, i, 
                                multiple_outputs,
                                separator,
                                step_data['sort_options'],
                                memcap),
                             callback=return_values.append
                        )
                pool.close()
                while len(return_values) != input_file_count:
                    try:
                        max_tuple = max(map(len, return_values))
                    except ValueError:
                        # return_values is empty
                        max_tuple = -1
                        pass
                    if max_tuple > 0:
                        # There are error tuples
                        errors = ['Streaming command "%s" failed: %s' % 
                                  (error[1], error[0]) for error 
                                  in return_values if len(error) == 2]
                        errors = [('%d) ' % (i + 1)) + error
                                    for i, error in enumerate(errors)]
                        iface.fail('\n'.join(errors),
                                   (job_flow[step_number:]
                                    if step_number != 0 else None))
                        failed = True
                        raise RuntimeError
                    iface.status('    Tasks completed: %d/%d'
                                 % (len(return_values), input_file_count))
                    time.sleep(0.2)
                try:
                    max_tuple = max(map(len, return_values))
                except ValueError:
                    # return_values is empty
                    max_tuple = -1
                    pass
                if max_tuple > 0:
                    # There are error tuples
                    errors = ['Streaming command "%s" failed: %s' % 
                              (error[1], error[0]) for error 
                              in return_values if len(error) == 2]
                    errors = [('%d) ' % (i + 1)) + error
                                for i, error in enumerate(errors)]
                    iface.fail('\n'.join(errors),
                               (job_flow[step_number:]
                                if step_number != 0 else None))
                    failed = True
                    raise RuntimeError
                iface.step('    Completed %s.'
                           % dp_iface.inflected(input_file_count, 'task'))
            # Really close open file handles in PyPy
            gc.collect()
            if not keep_intermediates:
                iface.status('    Deleting temporary files....')
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
        if not keep_last_output:
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
    except Exception, GeneratorExit:
        # GeneratorExit added just in case this happens on modifying code
        if not failed:
            time.sleep(0.2)
            if 'step_number' in locals():
                iface.fail(steps=(job_flow[step_number:]
                            if step_number != 0 else None))
            else:
                iface.fail()
        raise
    except KeyboardInterrupt, SystemExit:
        if 'step_number' in locals():
            iface.fail(steps=(job_flow[step_number:]
                        if step_number != 0 else None),
                        opener='*****Terminated*****')
        else:
            iface.fail()
        if 'pool' in locals():
            pool.terminate()
            pool.join()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, 
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
            '-m', '--memcap', type=int, required=False, default=(1024*300),
            help='Maximum amount of memory to use per UNIX sort instance.'
        )
    parser.add_argument(
            '-p', '--num-processes', type=int, required=False, default=1,
            help='Number of subprocesses to open at once.'
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
    dp_iface.add_args(parser)
    args = parser.parse_args(sys.argv[1:])
    run_simulation(args.branding, args.json_config, args.force,
                    args.memcap, args.num_processes, args.separator,
                    args.keep_intermediates, args.keep_last_output,
                    args.log)