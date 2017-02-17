#!/usr/bin/env python
"""
hadoop_runner.py
Part of Dooplicity framework

Outputs a bash script with Hadoop Streaming job flow.
FUNCTIONALITY IS IDIOSYNCRATIC;

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
import json

def add_args(parser):
    """ Adds args relevant to Hadoop runner.

        parser: object of type parser.ArgumentParser

        No return value.
    """
    parser.add_argument('-d', '--hadoop-path', type=str,
                        required=True, default=None, help=('Location \
                        of hadoop for running the hadoop job \
                        flow. Type `which hadoop` to find the \
                        location.'))
    parser.add_argument('-j', '--job-flow', required=True,
                        default=None, help=('File that contains Hadoop \
                        job flow described in json format.')  )
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
                        exit. ' 'Default is false.'))

def extract_steps_input_output(hadoop_path, job_flow):
    # TODO: documentation
    # Return a dictionary of the inputs, a set of outputs, and a list of hadoop commands.
    steps = []
    hadoop_step = hadoop_path
    input_last_seen = {}
    all_outputs = set()
    currently_input = False
    currently_output = False

    for step in job_flow["Steps"]:
        for arg in step["HadoopJarStep"]["Args"]:
            hadoop_step += " " + arg
            # If currently seeing an input, record the latest step it
            # was seen, which is the current step.
            if currently_input:
                for step_input in arg.split(','):
                    input_last_seen[step_input] = len(steps)
                currently_input = False
            # If currently seeing an output, add it to the set.
            elif currently_output:
                print arg
                all_outputs.add(arg)
                currently_output = False
            # If current argument is "-input", the next argument will
            # be an input directory.
            if arg == "-input":
                currently_input = True
            elif arg == "-output":
                currently_output = True
        steps.append(hadoop_step)
        hadoop_step = hadoop_path

    return input_last_seen, all_outputs, steps

def add_delete_intermediate_steps(steps, input_last_seen):
    # TODO: documentation
    # TODO: write function
    return steps

def get_hadoop_streaming_command(hadoop_path, job_flow, keep_intermediates):
    # TODO: documentation
    input_last_seen, all_outputs, steps = extract_steps_input_output(hadoop_path, job_flow)
    if not keep_intermediates:
        steps = add_delete_intermediate_steps(steps, input_last_seen)

    return steps

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    add_args(parser)
    args = parser.parse_args(sys.argv[1:])
    with open(args.job_flow, 'r') as job_flow_file:
        job_flow = json.load(job_flow_file)
    hadoop_commands = get_hadoop_streaming_command(
                             args.hadoop_path, job_flow, args.keep_intermediates)
    if args.print_command:
          print hadoop_commands
    if args.output_path:
        with open(args.output_path, 'w') as hadoop_command_file:
            for hadoop_command in hadoop_commands:
                hadoop_command_file.write(hadoop_command)
                hadoop_command_file.write("\n")
