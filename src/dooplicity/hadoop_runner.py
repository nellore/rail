#!/usr/bin/env python
"""
hadoop_runner.py
Part of Dooplicity framework

Outputs a bash script with Hadoop Streaming job flow.
FUNCTIONALITY IS IDIOSYNCRATIC;

Licensed under the MIT License:

Copyright (c) 2017 Jeena Lee and Abhi Nellore.

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
import json


class HadoopInputLastSeen(object):
    """ Stores the hadoop input file and its last seen step.

        Attributes:
            hadoop_input: A string of a hadoop input file for a hadoop step.
            step: An integer for a step the input file was last seen.
    """
    def __init__(self, hadoop_input, step):
        """ Inits HadoopInputLastSeen with input and step. """
        self.hadoop_input = hadoop_input
        self.step = step


def add_args(parser):
    """ Adds args relevant to Hadoop runner.

        parser: object of type parser.ArgumentParser

        No return value.
    """
    parser.add_argument('-d', '--hadoop-path', type=str, default=None,
                        required=True, help=('Location of hadoop for running \
                        the hadoop job flow. Type `which hadoop` to find the \
                        location.'))
    # TODO
    # When the json formatter is correctly written, hadoop
    # streaming path will be written in the json file. Therefore, `-s`
    # argument can be removed.
    parser.add_argument('-s', '--hadoop-streaming-path', type=str, default=None,
                        required=True, help=('Location of hadoop streaming jar \
                        for running the hadoop job flow.'))
    parser.add_argument('-j', '--job-flow', default=None, required=True,
                        help=('File that contains Hadoop job flow described in \
                        json format.'))
    parser.add_argument('-k', '--keep-intermediates',
                        action='store_const', const=True,
                        default=False, help='Keeps all intermediate \
                        output.')
    parser.add_argument('-o', '--output-path', type=str,
                        required=False, default=None, help=('Output \
                        path to save the hadoop streaming jar \
                        command.'))
    parser.add_argument('-p', '--print-command', action='store_const',
                        const=True, required=False, default=False,
                        help=('Print the hadoop command line and \
                        exit. ' 'Default is False.'))
    parser.add_argument('-t', '--skip-trash', action='store_const',
                        const=True, required=False, default=True,
                        help=('When deleting files, skip trash and \
                        delete permanantly. Default is True.'))
    parser.add_argument('--test', action='store_const', const=True,
                        default=False,
                        help='Run unit tests; DOES NOT NEED INPUT FROM STDIN, \
                        AND DOES NOT WRITE EXONS AND JUNCTIONS TO STDOUT')

def extract_steps_input_output(hadoop_path, hadoop_streaming_path, job_flow):
    """ Extracts hadoop steps, inputs, and outputs from a hadoop job flow.

        hadoop_path: Path to hadoop given by user.
        hadoop_streaming_path: Path to hadoop streaming jar given by user.
        job_flow: JSON object that describes a hadoop job flow.

        Return values:
            steps: A list of hadoop commands without
                   intermediate data deletion steps.
            input_last_seen: A list of HadoopInputLastSeen objects.
            all_outputs: A set of all outputs found in job_flow.
    """
    steps = []
    hadoop_step = hadoop_path + " jar "
    input_last_seen = []
    all_outputs = set()
    # `state` is the type of argument currently it's reading.
    state = None

    for step in job_flow["Steps"]:
        # TODO
        # In the future when dooplicity creates correct json
        # workflow with the appropriate jar file location, remove
        # manually giving hadoop_streaming_path as an argument to
        # hadoop_runner.py.
        if step["HadoopJarStep"]["Jar"] == "":
            hadoop_step += hadoop_streaming_path
        else:
            hadoop_step += step["HadoopJarStep"]["Jar"]

        # Process the arguments for a hadoop streaming step.
        for arg in step["HadoopJarStep"]["Args"]:
            # For `-mapper` and `-reducer` args, they must be wrapped
            # in `"`. Otherwise, append the arg as it is.
            if state == "-D":
                hadoop_step += '"' + arg + '"'
                state = None
            elif state == "-mapper" or state == "-reducer":
                hadoop_step += ' "' + arg + '"'
                state = None
            else:
                hadoop_step += ' ' + arg

            # If current state is input, record the latest step the
            # input was seen, which is the current step.
            if state == "-input":
                for step_input in arg.split(','):
                    input_last_seen = update_input_last_seen(
                        step_input, len(steps), input_last_seen)
                state = None
            # If current state is output, add output to the output set.
            elif state == "-output":
                all_outputs.add(arg)
                state = None

            # TODO
            # This will be removed when the json formatter is
            # improved. The inputformat `com.twitter.elephantbird. \
            # mapred.input.DeprecatedCombineLzoTextInputFormat`, when
            # run on hadoop, appends the line number in front of the
            # line, increasing the line length by 1. This disrupts the
            # downstream process, and therefore the first column of a
            # line should be removed at the map step if and only if
            # `-mapper` is `cat`. `-mapper` should be changed to `cut -f2-`.
            elif state == "-inputformat":
                if arg == "com.twitter.elephantbird.mapred.input.DeprecatedCombineLzoTextInputFormat":
                    if '-mapper "cat"' in hadoop_step:
                        hadoop_step = hadoop_step.replace('-mapper "cat"',
                                                          '-mapper "cut -f2-"')
                state = None

            # Depending on the current argument, set the `state`.
            if arg == "-input":
                state = "-input"
            elif arg == "-output":
                state = "-output"
            elif arg == "-D":
                state = "-D"
            elif arg == "-mapper":
                state = "-mapper"
            elif arg == "-reducer":
                state = "-reducer"
            elif arg == "-inputformat":
                state = "-inputformat"
        steps.append(hadoop_step)

        hadoop_step = hadoop_path + " jar "

    # `input_last_seen` is reversely sorted so that we can add steps
    # for deleting intermediate data later if needed.
    input_last_seen = sorted(input_last_seen,
                             key=lambda x: x.step, reverse=True)

    return steps, input_last_seen, all_outputs


def update_input_last_seen(step_input, step_length, input_last_seen):
    """ Updates input_last_seen to correctly reflect its last seen step.

        step_input: An input for a step.
        step_length: Total number of steps that were seen so far.
        input_last_seen: A list of HadoopInputLastSeen objects.

        Return values:
            input_last_seen: A list of HadoopInputLastSeen objects.
    """
    seen_before = False
    for input_step_pair in input_last_seen:
        if input_step_pair.hadoop_input == step_input:
            input_step_pair.step = step_length
            seen_before = True
    if not seen_before:
        input_last_seen.append(HadoopInputLastSeen(step_input, step_length))

    return input_last_seen


def add_delete_intermediate_steps(hadoop_path, steps,
                                  input_last_seen, all_outputs, skip_trash):
    """ Adds steps to delete intermediates to a list of hadoop commands.

        hadoop_path: Path to hadoop given by user.
        steps: A list of hadoop commands without
              intermediate data deletion steps.
        input_last_seen: A list of HadoopInputLastSeen objects.
        all_outputs: A set of all outputs found in job_flow.
        skip_trash: A boolean for whether to skip trash
                   when deleting intermediates from HDFS.

        Return values:
            steps: A list of hadoop commands with
                  intermediate data deletion steps.
    """
    delete_command = hadoop_path + " fs -rm -r "
    if skip_trash:
        delete_command += "-skipTrash "

    for to_remove in input_last_seen:
        if to_remove.hadoop_input not in all_outputs:
            '''Remove directory only if it's an -output of some
            step and an -input of another step.'''
            continue
        else:
            delete_command += to_remove.hadoop_input
            steps.insert(to_remove.step + 1, delete_command)
        delete_command = hadoop_path + " fs -rm -r "
        if skip_trash:
            delete_command += "-skipTrash "

    return steps


def get_hadoop_streaming_command(hadoop_path,
                                 hadoop_streaming_path,
                                 job_flow,
                                 keep_intermediates,
                                 skip_trash):
    """ Add steps to delete intermediates to a list of hadoop commands.

        hadoop_path: Path to hadoop given by user.
        job_flow: JSON object that describes a hadoop job flow.
        keep_intermediates: A boolean for whether to store
                            intermediate data.
        skip_trash: A boolean for whether to skip trash
                   when deleting intermediates from HDFS.

        Return values:
            steps: A list of hadoop commands with or without
                  intermediate data deletion steps depending on
                  keep_intermediates.
    """
    steps, input_last_seen, all_outputs = extract_steps_input_output(
        hadoop_path, hadoop_streaming_path, job_flow)
    if not keep_intermediates:
        steps = add_delete_intermediate_steps(hadoop_path,
                                              steps,
                                              input_last_seen,
                                              all_outputs,
                                              skip_trash)

    return steps

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    add_args(parser)
    args = parser.parse_args(sys.argv[1:])

if __name__ == '__main__' and not args.test:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    with open(args.job_flow) as job_flow_file:
        job_flow = json.load(job_flow_file)
    hadoop_commands = get_hadoop_streaming_command(args.hadoop_path,
                                                   args.hadoop_streaming_path,
                                                   job_flow,
                                                   args.keep_intermediates,
                                                   args.skip_trash)
    if args.print_command:
        print hadoop_commands
    if args.output_path:
        with open(args.output_path, 'w') as hadoop_command_file:
            print >>hadoop_command_file, (
"""#!/usr/bin/env bash

set -e"""
                )
            for hadoop_command in hadoop_commands:
                print >>hadoop_command_file, hadoop_command
elif __name__ == '__main__':
  # Test units
    del sys.argv[1:] # Don't choke on extra command-line parameters
    import unittest
    import tempfile
    import shutil

    class TestHadoopCommandOutput(unittest.TestCase):
        """ Unit tests for getting hadoop streaming command from a json job flow. """
        def setUp(self):
            self.job_flow = {'Steps':
                         [{'HadoopJarStep':
                           {'Args':
                            ['-D', 'mapreduce.job.reduces=0',
                              '-input', 'file_0', '-output', 'file_1'], 'Jar': ''},
                           'Name': 'First step'
                         },
                          {'HadoopJarStep':
                           {'Args':
                            ['-D', 'mapreduce.job.reduces=1', '-D',
                             'stream.num.map.output.key.fields=1',
                             '-input', 'file_1', '-output', 'file_2'], 'Jar': ''},
                           'Name': 'Second step'
                          },
                          {'HadoopJarStep':
                           {'Args':
                            ['-D', 'mapreduce.job.reduces=2',
                             '-input', 'file_1', '-output', 'file_3',
                             '-mapper', 'cat'], 'Jar': ''},
                           'Name': 'Third step'
                              },
                          {'HadoopJarStep':
                           {'Args':
                            ['-D', 'mapreduce.job.reduces=3',
                             '-input', 'file_3,file_2', '-output', 'file_4'], 'Jar': ''},
                           'Name': 'Fourth step'
                          }
                         ]}

            self.job_flow_map = {'Steps':
                             [{'HadoopJarStep':
                               {'Args':
                                ['-D', 'mapreduce.job.reduces=0',
                                 '-input', 'file_0', '-output', 'file_1',
                                 '-mapper', 'cat', '-inputformat',
                                 'com.twitter.elephantbird.mapred.input.DeprecatedCombineLzoTextInputFormat'
                            ], 'Jar': ''},
                               'Name': 'First step'
                             },
                              {'HadoopJarStep':
                               {'Args':
                                ['-D', 'mapreduce.job.reduces=1', '-D',
                                 'stream.num.map.output.key.fields=1',
                                 '-input', 'file_1', '-output', 'file_2',
                                 '-mapper', 'count.py', '-inputformat',
                                 'com.twitter.elephantbird.mapred.input.DeprecatedCombineLzoTextInputFormat'
                                ], 'Jar': ''},
                               'Name': 'Second step'
                              }
                             ]}

        def tearDown(self):
            self.job_flow = None
            self.job_flow_map = None

        def test_keep_intermediate_and_skip_trash(self):
            hadoop_commands = get_hadoop_streaming_command("hadoop",
                                                       "hadoop-streaming-jar",
                                                       self.job_flow,
                                                       keep_intermediates=False,
                                                       skip_trash=True)
            expected_hadoop_commands = [
            'hadoop jar hadoop-streaming-jar -D"mapreduce.job.reduces=0" \
-input file_0 -output file_1', \
            'hadoop jar hadoop-streaming-jar -D"mapreduce.job.reduces=1" \
-D"stream.num.map.output.key.fields=1" -input file_1 -output file_2', \
            'hadoop jar hadoop-streaming-jar -D"mapreduce.job.reduces=2" \
-input file_1 -output file_3 -mapper "cat"', \
            'hadoop fs -rm -r -skipTrash file_1', \
            'hadoop jar hadoop-streaming-jar -D"mapreduce.job.reduces=3" \
-input file_3,file_2 -output file_4', \
            'hadoop fs -rm -r -skipTrash file_2', \
            'hadoop fs -rm -r -skipTrash file_3']
            self.assertEqual(hadoop_commands, expected_hadoop_commands)

        def test_keep_intermediate_and_not_skip_trash(self):
            hadoop_commands = get_hadoop_streaming_command("hadoop",
                                                       "hadoop-streaming-jar",
                                                       self.job_flow,
                                                       keep_intermediates=False,
                                                       skip_trash=False)
            expected_hadoop_commands = [
            'hadoop jar hadoop-streaming-jar -D"mapreduce.job.reduces=0" \
-input file_0 -output file_1', \
            'hadoop jar hadoop-streaming-jar -D"mapreduce.job.reduces=1" \
-D"stream.num.map.output.key.fields=1" -input file_1 -output file_2', \
            'hadoop jar hadoop-streaming-jar -D"mapreduce.job.reduces=2" \
-input file_1 -output file_3 -mapper "cat"', \
            'hadoop fs -rm -r file_1', \
            'hadoop jar hadoop-streaming-jar -D"mapreduce.job.reduces=3" \
-input file_3,file_2 -output file_4', \
            'hadoop fs -rm -r file_2', \
            'hadoop fs -rm -r file_3']
        self.assertEqual(hadoop_commands, expected_hadoop_commands)

        def test_basic_hadoop_command_output(self):
            hadoop_commands = get_hadoop_streaming_command("hadoop",
                                                       "hadoop-streaming-jar",
                                                       self.job_flow,
                                                       keep_intermediates=True,
                                                       skip_trash=False)
            expected_hadoop_commands = [
            'hadoop jar hadoop-streaming-jar -D"mapreduce.job.reduces=0" \
-input file_0 -output file_1', \
            'hadoop jar hadoop-streaming-jar -D"mapreduce.job.reduces=1" \
-D"stream.num.map.output.key.fields=1" -input file_1 -output file_2', \
            'hadoop jar hadoop-streaming-jar -D"mapreduce.job.reduces=2" \
-input file_1 -output file_3 -mapper "cat"', \
            'hadoop jar hadoop-streaming-jar -D"mapreduce.job.reduces=3" \
-input file_3,file_2 -output file_4']
        self.assertEqual(hadoop_commands, expected_hadoop_commands)

        def test_mapper_hadoop_command_output(self):
            hadoop_commands = get_hadoop_streaming_command("hadoop",
                                                       "hadoop-streaming-jar",
                                                       self.job_flow_map,
                                                       keep_intermediates=True,
                                                       skip_trash=False)
            expected_hadoop_commands = [
            'hadoop jar hadoop-streaming-jar -D"mapreduce.job.reduces=0" \
-input file_0 -output file_1 -mapper "cut -f2-" -inputformat \
com.twitter.elephantbird.mapred.input.DeprecatedCombineLzoTextInputFormat',
            'hadoop jar hadoop-streaming-jar -D"mapreduce.job.reduces=1" \
-D"stream.num.map.output.key.fields=1" -input file_1 -output file_2 \
-mapper "count.py" -inputformat \
com.twitter.elephantbird.mapred.input.DeprecatedCombineLzoTextInputFormat']
            self.assertEqual(hadoop_commands, expected_hadoop_commands)

    unittest.main()