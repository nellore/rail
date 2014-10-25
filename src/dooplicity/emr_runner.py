#!/usr/bin/env python
"""
emr_runner.py
Part of Dooplicity framework

Runs JSON-encoded Hadoop Streaming job flow on Elastic MapReduce. Interface
should eventually match that of emr_simulator.py as closely as possible.
Uses AWS CLI to talk to S3.

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
import interface as dp_iface
import sys
import json
from collections import OrderedDict, defaultdict
import argparse
import webbrowser
import ansibles as ab

_aws_regions = set(['us-east-1', 'us-west-2', 'us-west-1', 'eu-west-1',
                    'ap-southeast-1', 'ap-southeast-2', 'ap-northeast-1',
                    'sa-east-1', 'us-gov-west-1'])
_emr_url = ('https://console.aws.amazon.com/elasticmapreduce/?region={region}'
            '#cluster-details:{job_flow_id}')

def add_args(parser):
    """ Adds args relevant to EMR job submission.

        parser: object of type parser.ArgumentParser

        No return value.
    """
    parser.add_argument(
            '-x', '--aws-exe', type=str, required=False, default='aws',
            help='Path to AWS CLI executable. Specifying this should be '
                 'unnecessary on all platforms if the AWS CLI is installed '
                 'properly, but the user may want to select a different '
                 'executable if multiple versions are installed on her '
                 'system.'
        )
    parser.add_argument(
            '-r', '--region', type=str, required=False, default=None,
            help='Amazon data center in which to run job flow. Defaults to ' \
                 'region from AWS CLI; specifying the region here overrides ' \
                 'any region in environment variables and from "--profile".'
        )
    parser.add_argument(
            '-i', '--profile', type=str, required=False, default='default',
            help='Name of AWS CLI profile from which to grab AWS access key ' \
                 'ID and secret access key. Defaults to "default" profile.'
        )
    parser.add_argument('-n', '--no-browser', action='store_const', const=True,
            default=False,
            help='Do not automatically open browser after job flow ' \
                 'submission.')

def run_job_flow(branding, json_config, force, no_browser=False,
                    aws_exe='aws', profile='default', region=None):
    """ Submits job flow to EMR. Should eventually monitor it.

        branding: text file with branding to print to screen before running
            job. This is where the name of a software package or ASCII art 
            can go.
        json_config: JSON configuration file. Google Getting Started with
            Amazon Elastic MapReduce for formatting information.
        force: True iff all existing directories should be erased when
            writing intermediates.
        no_browser: True iff webpage with status of job shouldn't automatically
            be opened in default browser.
        aws_exe: path to AWS CLI executable
        profile: name of AWS CLI profile from which to grab AWS access key ID
            and secret access key.
        region: Amazon data center in which to run job flow. If None, uses
            region from AWS CLI.
        No return value.
    """
    iface = dp_iface.DooplicityInterface(
                        branding=branding,
                        opener='Started job flow submission script on {time}.'
                    )
    try:
        '''Should ultimately check to see which regions are available to user
        w/ API.'''
        aws_ansible = ab.AWSAnsible(profile=profile, region=region,
                                    valid_regions=_aws_regions)
        s3_ansible = ab.S3Ansible(aws_exe=aws_exe, profile=profile,
                                    iface=iface)
        # Serialize JSON configuration
        if json_config is not None:
            with open(json_config) as json_stream:
                full_payload = json.load(json_stream)
        else:
            full_payload = json.load(sys.stdin)
        try:
            job_flow = full_payload['Steps']
        except Exception:
            raise RuntimeError(
                    'Input JSON not in proper format. Ensure that the JSON '
                    'object has a Steps key.'
                )
        step_count = len(job_flow)
        steps = OrderedDict()
        # Check steps for requred data
        required_data = set(['input', 'output', 'mapper', 'reducer'])
        bad_output_data = []
        try:
            for step in job_flow:
                if 'hadoop-streaming' \
                    not in step['HadoopJarStep']['Jar'].lower():
                    '''Don't check command-line parameters if not a 
                    Hadoop Streaming step.'''
                    continue
                step_args = {}
                for j in xrange(0, len(step['HadoopJarStep']['Args']), 2):
                    arg_name = step['HadoopJarStep']['Args'][j][1:].strip()
                    if arg_name == 'input':
                        try:
                            step_args['input'] = ','.join(
                                    [step['HadoopJarStep']['Args'][j+1],
                                     step_args['input']]
                                )
                        except KeyError:
                            step_args['input'] \
                                = step['HadoopJarStep']['Args'][j+1]
                    elif arg_name in required_data:
                        step_args[step['HadoopJarStep']['Args'][j][1:]] \
                            = step['HadoopJarStep']['Args'][j+1].strip()
                steps[step['Name']] = step_args
        except (KeyError, IndexError):
            iface.fail(
                    'JSON file not in proper format. Ensure '
                    'that each StepConfig object has a HadoopJarStep '
                    'object with an Args array and a Name string.'
                )
            raise
        iface.step('Read job flow from input JSON.')
        iface.status('Checking that output directories on S3 are writable...')
        '''Check steps for required Hadoop Streaming command-line parameters
        and for whether outputs are writable.'''
        missing_data = defaultdict(list)
        identity_steps = []
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
                    and ab.Url(step_data['output']).is_s3 \
                         and s3_ansible.is_dir(step_data['output']):
                    bad_output_data.append(step)
            if 'mapper' in steps[step] and 'reducer' in step_data \
                and step_data['mapper'] in identity_mappers \
                and step_data['reducer'] in identity_reducers:
                identity_steps.append(step)
        errors = []
        if missing_data:
            errors.extend(['Step "%s" is missing required parameter(s) "%s".' % 
                                (step, ', '.join(missing_data[step]))
                                for step in missing_data])
        if bad_output_data:
            errors.extend(['Output "directory" "%s" of step "%s" already '
                           'exists, and --force was not invoked to permit '
                           'overwriting it.' 
                           % (steps[step]['output'], step)
                           for step in bad_output_data])
        if identity_steps:
            errors.extend([('Step "%s" has both an identity mapper and an '
                            'identity reducer, which is redundant. Remove the '
                            'step before trying again.') % step])
        if errors:
            raise RuntimeError('\n'.join([('%d) ' % (i+1)) + error
                                    for i, error in enumerate(errors)]))
        iface.step('Verified that output directories on S3 are writable.')
        '''Remove intermediate directories and create buckets if they
        don't already exist so EMR doesn't choke when writing output. The
        idea here is to make a bucket feel like an ordinary directory to
        the user.'''
        buckets = set()
        for step in steps:
            step_data = steps[step]
            if ab.Url(step_data['output']).is_s3:
                s3_ansible.remove_dir(steps[step]['output'])
                buckets.add(ab.bucket_from_url(step_data['output']))
        for bucket in buckets:
            try:
                s3_ansible.create_bucket(bucket)
            except Exception as e:
                raise RuntimeError(('Bucket %s already exists on S3. Change '
                                    'affected output directories in job flow '
                                    'and try again. The more distinctive the '
                                    'name chosen, the less likely it\'s '
                                    'already a bucket on S3. It may be '
                                    'easier to create a bucket first using '
                                    'the web interface and use its name + a '
                                    'subdirectory as the output directory.')
                                    % bucket)
        iface.step('Set up output directories on S3.')
        iface.status('Submitting job flow...')
        job_flow_response = aws_ansible.post_request(full_payload)
        json_response = json.load(job_flow_response)
        if 'JobFlowId' not in json_response:
            raise RuntimeError('Job flow submission failed. Server returned '
                               'the following response: %s.'
                                % json.dumps(json_response, sort_keys=True,
                                    indent=4, separators=(',', ': ')))
        else:
            job_flow_id = json_response['JobFlowId']
            iface.step(('Submitted job flow.\n'
                        '*****Job flow ID is %s .*****\n' % job_flow_id) + 
                        ('*****Submission can be monitored at %s .*****' %
                            _emr_url.format(region=aws_ansible.region,
                                            job_flow_id=job_flow_id)))
            if not no_browser:
                webbrowser.open(_emr_url.format(region=aws_ansible.region,
                                                job_flow_id=job_flow_id),
                                    new=2) # Open in new tab
                iface.step('Opening URL in default browser.')
        iface.done(closer=('Finished job flow submission script '
                           'on {date}. Run time was {length} seconds.')
                )
    except (Exception, GeneratorExit):
         # GeneratorExit added just in case this happens on modifying code
        iface.fail()
        raise
    except (KeyboardInterrupt, SystemExit):
        if 'job_flow_id' not in locals():
            iface.fail(opener='*****Terminated without submitting job.*****',
                        middler=('End time was {date}. '
                                 'Script was running for {length} seconds.'))
        else:
            iface.fail(opener='*****Terminated after submitting job.*****',
                        middler=('End time was {date}. '
                                 'Script was running for {length} seconds.'))
        raise

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, 
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    add_args(parser)
    dp_iface.add_args(parser)
    args = parser.parse_args(sys.argv[1:])
    run_job_flow(args.branding, args.json_config, args.force, args.no_browser,
                    args.aws_exe, args.profile, args.region)
