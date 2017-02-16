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
    parser.add_argument(
            '-d', '--hadoop-path', type=str, required=True, default=None,
            help=('Location of hadoop for running the hadoop job flow. '
                  'Type `which hadoop` to find the location.')
    )
    parser.add_argument(
            '-j', '--job-flow', required=True, default=None,
            help=('File that contains Hadoop job flow described in json format.')
    )
    parser.add_argument(
        '-p', '--print-command',
        action='store_const', const=True, required=False,
        default=False, help=('Print the hadoop command line and \
        exit. ' 'Default is false.')
    )
    parser.add_argument(
            '-o', '--output-path', type=str, required=False, default=None,
            help=('Output path to save the hadoop streaming jar command.')
    )

# Scan the entire job flow
# Create a dictionary for each input file. Value will be the latest step it was seen
# Same time, create hadoop commands and add to a list
# Later, create deletion step and insert into the list of commands

def extract_hadoop_commands_input_dir(hadoop_path, job_flow):
    # TODO doc
    # Return a dictionary of the input files and a list of hadoop commands.
    hadoop_commands = []
    hadoop_command = hadoop_path
    input_dict = {}
    currently_input = False

    for step in job_flow["Steps"]:
        for arg in step["HadoopJarStep"]["Args"]:
            hadoop_command += " " + arg
            # If currently seeing input, record the latest step it was seen.
            if currently_input:
                input_dict[arg] = len(hadoop_commands)
                currently_input = False
            # If current argument is "-input", the next argument will be an input directory.
            if arg == "-input":
                currently_input = True
        hadoop_commands.append(hadoop_command)
        hadoop_command = hadoop_path

    return input_dict, hadoop_commands

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    add_args(parser)
    args = parser.parse_args(sys.argv[1:])
    # TODO:
    # Handle -p argument
    # Handle -o argument
    with open(args.job_flow, 'r') as job_flow_file:
        job_flow = json.load(job_flow_file)
    job_flow_file.close()
    dictionary, hadoop_commands = extract_hadoop_commands_input_dir(
                                      args.hadoop_path, job_flow)

    if args.print_command:
          print hadoop_commands
