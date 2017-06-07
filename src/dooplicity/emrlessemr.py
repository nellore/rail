#!/usr/bin/env python
"""
emrlessemr.py

Writes shells scripts for a) knitting a set of EC2 instances, each running an
EMR AMI, into a Hadoop cluster without EMR; b) running bootstraps from EMR JSON
on each instance; and c) running a job flow from EMR JSON. All scripts should 
be run from what will be the master instance.
"""
import argparse
import re
import json
import os
import errno

if __name__ == '__main__':
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-j', '--json', type=str, required=False,
            default='-',
            help='path to job flow JSON or - if read from stdin'
        )
    parser.add_argument('-i', '--ips', type=str, required=False,
            default='-',
            help='path to EC2 IP list for core intances or - if read from '
                 'stdin'
        )
    parser.add_argument('-o', '--output', type=str, required=False,
            default='.'
            help='directory in which to write shell scripts'
        )
    parser.add_argument('--aws', type=str, required=False,
            default='aws'
            help='path to AWS CLI executable; must be at same path on every '
                 'EC2 instance'
        )
    import re
    if not sys.stdin.isatty():
        # Try to grab job flow and IPs from stdin
        tokens = re.split('({.+})', sys.stdin.read())
        ips = set([ip.strip() for ip
                        in tokens[0].splitlines() + tokens[-1].splitlines()])
        try:
            parsed_json = json.loads(tokens[1])
        except IndexError:
            parsed_json = None
    # Read JSON and IPs from files
    if args.json == '-':
        if parsed_json == None:
            raise RuntimeError(
                    'Job flow JSON read from stdin was not properly '
                    'formatted; make sure it begins and ends with braces '
                    'before trying again.'
                )
    else:
        with open(args.json) as json_stream:
            parsed_json = json.load(json_stream)
    if args.ips == '-':
        if not args.ips:
            raise RuntimeError(
                    'No EC2 IPs were found in stdin.'
                )
    else:
        with open(args.ips) as ip_stream:
            ips = set(ip_stream.read().splitlines())
    
    args.output = os.path.abspath(os.path.expandvars(args.output))
    try:
        os.path.makedirs(args.output)
    except OSError as e:
        # Forgive existing directories
        if e.errno != errno.EEXIST:
            raise

    os.chdir(args.output)
    with open('start_bootstraps.sh', 'w') as script_stream:
        parsed_json['BootstrapActions']
