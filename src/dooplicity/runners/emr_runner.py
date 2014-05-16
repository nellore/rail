#!/usr/bin/env python
"""
emr_runner.py
Part of Dooplicity framework

Runs JSON-encoded Hadoop Streaming job flow on Elastic MapReduce. Interface
should eventually match that of hadoop_simulator.py as closely as possible.
Uses s3cmd to talk to S3.

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
import dooplicity_interface as dp_iface
import subprocess
import sys
import os
import time
# For signing requests as per AWS API
import hashlib
import hmac
import binascii
import urllib
import urllib2
import json
from collections import OrderedDict, defaultdict
from dooplicity_version import version as _version
import argparse
import webbrowser
import site
import math
from functools import wraps

site.addsitedir(os.path.dirname(os.path.dirname(__file__)))
import dooplicity as dp

_aws_regions = set(['us-east-1', 'us-west-2', 'us-west-1', 'eu-west-1',
                    'ap-southeast-1', 'ap-southeast-2', 'ap-northeast-1',
                    'sa-east-1', 'us-gov-west-1'])
_emr_url = ('https://console.aws.amazon.com/elasticmapreduce/'
            '#cluster-details:{job_flow_id}')

def clean_url(url):
    """ Tacks an s3:// onto the beginning of a URL if necessary. 

        url: an S3 URL, or what should be one

        Return value: S3 URL
    """
    if url[:6] == 's3n://':
        return ('s3://' + url[6:])
    elif url[:5] != 's3://':
        return ('s3://' + url)
    else:
        return url

def bucket_from_url(url):
    """ Extracts a bucket name from an S3 URL.

        url: an S3 URL, or what should be one

        Return value: bucket name from S3 URL
    """
    if url[:6] == 's3n://':
        start_index = 6
    elif url[:5] == 's3://':
        start_index = 5
    elif url[0] == '/':
        start_index = 1
    else:
        start_index = 0
    while url[start_index] == '/':
        start_index += 1
    return url[start_index:start_index+url[start_index:].index('/')]

class S3Ansible:
    """ Permits simple interactions with S3 via s3cmd. 

        Another option is to use Boto, but it comes with a lot of unused
        classes --- tools never touched by Dooplicity/Rail. Further, my
        preference is to integrate tools with good command-line APIs so
        languages can be mixed and matched. -AN
    """

    def __init__(self, s3cmd_exe='s3cmd', config=None):
        """ Initializes S3Ansible.

            s3cmd_exe: path to s3cmd executable
            config: s3cmd config file, if specified by user
        """
        if config is None:
            extra_args = ''
        else:
            extra_args = '-c %s' % config
        self.s3cmd = ' '.join([s3cmd_exe, extra_args])
        self._osdevnull = open(os.devnull, 'w')

    def is_bucket(self, bucket):
        """ Checks whether a bucket exists.

            bucket: string with bucket name to check; can have an s3:// at
                the beginning or not

            Return value: True iff bucket exists; else False.
        """
        s3cmd_command = ' '.join([self.s3cmd, 'info', clean_url(bucket)])
        s3cmd_process = subprocess.Popen(
                            s3cmd_command,
                            shell=True,
                            bufsize=-1,
                            stderr=subprocess.PIPE,
                            stdout=self._osdevnull
                        )
        s3cmd_error = s3cmd_process.stderr.read().strip()
        s3cmd_return = s3cmd_process.wait()
        if s3cmd_return != 0:
            raise RuntimeError('s3cmd failed with exit level %d'
                                % s3cmd_return)
        if s3cmd_error:
            if 'does not exist' in s3cmd_error:
                return False
            else:
                raise RuntimeError('s3cmd says, "%s".' % s3cmd_error)
        return True

    def create_bucket(self, bucket):
        """ Creates a new bucket only if it doesn't already exist.

            bucket: string with bucket name to create; can have an s3:// at
                the beginning or not.

            No return value.
        """
        if not self.is_bucket(bucket):
            s3cmd_command = ' '.join([self.s3cmd, 'mb', clean_url(bucket)])
            s3cmd_process = subprocess.Popen(
                                s3cmd_command,
                                shell=True,
                                bufsize=-1,
                                stderr=subprocess.PIPE,
                                stdout=self._osdevnull
                            )
            s3cmd_error = s3cmd_process.stderr.read().strip()
            s3cmd_return = s3cmd_process.wait()
            if s3cmd_return != 0:
                raise RuntimeError('s3cmd failed with exit level %d' 
                                    % s3cmd_return)
            if s3cmd_error:
                raise RuntimeError('s3cmd says, "%s".' % s3cmd_error)
        '''Can check for the word 'created' in s3cmd's stdout output, but this
        may change from version to version of s3cmd.'''

    def is_dir(self, path):
        """ Checks whether a directory on S3 exists.

            path: directory name on S3, including bucket

            Return value: True iff directory exists; else False
        """
        cleaned_url = clean_url(path)
        while cleaned_url[-1] == '/':
            cleaned_url = cleaned_url[:-1]
        s3cmd_command = ' '.join([self.s3cmd, 'ls', cleaned_url])
        s3cmd_process = subprocess.Popen(
                            s3cmd_command,
                            shell=True,
                            bufsize=-1,
                            stderr=subprocess.PIPE,
                            stdout=subprocess.PIPE
                        )
        s3cmd_error = s3cmd_process.stderr.read().strip()
        to_return = False
        for line in s3cmd_process.stdout:
            line = line.strip()
            if line[:3] == 'DIR' and cleaned_url in line:
                to_return = True
        s3cmd_return = s3cmd_process.wait()
        if to_return:
            return True
        if s3cmd_return != 0:
            raise RuntimeError('s3cmd failed with exit level %d' 
                                % s3cmd_return)
        if s3cmd_error:
            if 'does not exist' in s3cmd_error:
                return False
            else:
                raise RuntimeError('s3cmd says, "%s".' % s3cmd_error)
        return False # Empty

    def remove_dir(self, path):
        """ Removes a directory on S3.

            path: directory name on S3, including bucket

            No return value.
        """
        cleaned_url = clean_url(path)
        if len(cleaned_url) < 6:
            raise RuntimeError('s3cmd says, "%s".' % s3cmd_error)
        while cleaned_url[-1] == '/':
            cleaned_url = cleaned_url[:-1]
        cleaned_url += '/'
        if self.is_dir(path):
            s3cmd_command = ' '.join([self.s3cmd, 'del --recursive',
                                        cleaned_url])
            s3cmd_process = subprocess.Popen(
                                s3cmd_command,
                                shell=True,
                                bufsize=-1,
                                stderr=subprocess.PIPE,
                                stdout=self._osdevnull
                            )
            s3cmd_error = s3cmd_process.stderr.read().strip()
            s3cmd_return = s3cmd_process.wait()
            if s3cmd_return != 0:
                raise RuntimeError('s3cmd failed with exit level %d' 
                                    % s3cmd_return)
            if s3cmd_error:
                raise RuntimeError('s3cmd says, "%s".' % s3cmd_error)

def aws_params_from_json(json_object, prefix=''):
    """ Parses JSON object to generate AWS query string params.

        If recursion depth is exceeded, there is no justice in the universe.
        This function may be useful for generating simple GET requests
        from JSON but is not currently in use in the code.

        json_object: subset of JSON object to parse.

        Return value: parsed URL dictionary.
    """
    params = {}
    if isinstance(json_object, dict):
        for key in json_object:
            for subkey, value in aws_params_from_json(
                                                        json_object[key],
                                                        prefix + key + '.'
                                                ).items():
                params[subkey] = value
    elif isinstance(json_object, list):
        for i, item in enumerate(json_object):
            for key, value in aws_params_from_json(
                                                    json_object[i],
                                                    prefix + ('member.%d.'
                                                                % (i+1))
                                                ).items():
                params[key] = value
    elif isinstance(json_object, basestring):
        params = {prefix[:-1] : json_object}
    return params

def keys_from_config_file(s3cmd_config_file):
    """ Extracts AWS key and AWS secret key from s3cmd config file.

        s3cmd_config_file: a file in the format of .s3cfg. See s3cmd
            documentation for more information.

        Return value: tuple (AWS key, AWS secret key)
    """
    access_key, secret_key = None, None
    with open(s3cmd_config_file) as s3cmd_config_stream:
        for line in s3cmd_config_stream:
            tokens = [token.strip() for token in line.split('=')]
            if tokens[0] == 'access_key':
                access_key = tokens[1]
            elif tokens[0] == 'secret_key':
                secret_key = tokens[1]
            if None not in [access_key, secret_key]:
                return access_key, secret_key
    raise RuntimeError(('%s is an invalid s3cmd configuration file. Make sure '
                        'the configuration file has the following two lines:\n'
                        'access_key = <your AWS access key ID>\n'
                        'secret_key = <your AWS secret access key>\n'
                        'Note that if s3cmd has been installed properly, no '
                        'configuration file needs to be specified; '
                        'Dooplicity uses the access keys entered on setup of '
                        's3cmd by default.') % s3cmd_config_file)

def sign(key, message):
    """ Helper function for outputting AWS signatures.

        See docstring of aws_signature() for more information.
    """
    return hmac.new(key, message.encode('utf-8'), hashlib.sha256).digest()

def aws_signature(string_to_sign, secret_key,
                    datestamp=time.strftime('%Y%m%d', time.gmtime()),
                    region='us-east-1', service='elasticmapreduce'):
    """ Outputs AWS Signature v4 for a string to sign from a secret key.

        Taken from 
        http://docs.aws.amazon.com/general/latest/gr/signature-v4-examples.html
        .

        string_to_sign: string to sign.
        secret_key: AWS secret key
        datestamp: time.gmtime(), a struct_time object

        Return value: precisely the signature; no need to URL encode it.
    """
    return binascii.hexlify(
        sign(
            sign(
                sign(
                    sign(
                        sign(
                            ('AWS4' + secret_key).encode('utf-8'),
                            datestamp),
                        region),
                    service),
                'aws4_request'),
                string_to_sign
            )
        )

def retry(ExceptionToCheck, tries=4, delay=3, backoff=2, logger=None):
    """ Retry calling the decorated function using an exponential backoff.

        http://www.saltycrane.com/blog/2009/11/
            trying-out-retry-decorator-python/
        Original from: http://wiki.python.org/moin/PythonDecoratorLibrary#Retry

        ExceptionToCheck: the exception to check. May be a tuple of
            exceptions to check
        tries: number of times to try (not retry) before giving up
        delay: initial delay between retries in seconds
        backoff: backoff multiplier e.g. value of 2 will double the delay
            each retry
        logger: logger to use. If None, print

        Return value: retry wrapper.
    """
    def deco_retry(f):
        @wraps(f)
        def f_retry(*args, **kwargs):
            mtries, mdelay = tries, delay
            while mtries > 1:
                try:
                    return f(*args, **kwargs)
                except ExceptionToCheck, e:
                    msg = '%s, Retrying in %d seconds...' % (str(e), mdelay)
                    if logger is not None:
                        logger.warning(msg)
                    time.sleep(mdelay)
                    mtries -= 1
                    mdelay *= backoff
            return f(*args, **kwargs)
        return f_retry  # true decorator
    return deco_retry

@retry(urllib2.URLError, tries=4, delay=3, backoff=2)
def urlopen_with_retry(request):
    """ Facilitates retries when opening URLs.

        request: urllib2.Request object

        Return value: file-like object of type urllib2.urlopen with
            response
    """
    return urllib2.urlopen(request)

class AWSAnsible:
    """ Interacts with AWS via GET requests. See AWS API documentation. 

        Uses Signature Version 4, the latest version as of 5/15/14.

        http://stackoverflow.com/questions/1088715/
        how-to-sign-amazon-web-service-requests-from-the-python-app-engine

        helped make this class.
    """
    def __init__(self, aws_access_key, aws_secret_access_key,
                    region='us-east-1',
                    base_uri='https://elasticmapreduce.amazonaws.com'):
        self._aws_access_key = aws_access_key
        self._aws_secret_access_key = aws_secret_access_key
        if base_uri[:8] == 'https://':
            base_suffix = base_uri[8:]
            self.base_prefix = base_uri[:8]
        elif base_uri[:7] == 'http://':
            base_suffix = base_uri[7:]
            self.base_prefix = base_uri[:7]
        else:
            raise RuntimeError('Base URI for GET request must begin with '
                               'http:// or https://.')
        split_base_suffix = [segment for segment in base_suffix.split('/')
                                if segment]
        self.service = split_base_suffix[0].partition('.')[0]
        self.host_header = '.'.join([region, split_base_suffix[0]])
        try:
            self.canonical_uri = '/' + '/'.join([urllib.quote_plus(segment) 
                for segment in split_base_suffix[1:]])
        except IndexError:
            # No segments
            self.canonical_uri = '/'
        self.region = region
        self.algorithm = 'AWS4-HMAC-SHA256'

    def post_request(self, json_config, target='ElasticMapReduce.RunJobFlow'):
        """ Submits JSON payload to AWS using Signature Version 4.

            Reference: http://docs.aws.amazon.com/general/latest/gr/
                sigv4-create-canonical-request.html

            json_config: file containing json payload
            target: the X-Amz-Target from the request header. Google for more
                information; for submitting job flows, this is
                ElasticMapReduce.RunJobFlow.

            Return value: file-like object of type urllib2.urlopen with
                response
        """
        with open(json_config) as json_stream:
            payload \
                = ''.join([line.strip() for line in json_stream.readlines()])
        hashed_payload = hashlib.sha256(payload).hexdigest()
        headers = {
                'Content-Type' : 'application/x-amz-json-1.1',
                'X-Amz-Target' : target,
                'Content-Length' : str(len(payload)),
                'User-Agent' : 'Dooplicity v%s' % _version,
                'Host' : self.host_header,
                'X-Amz-Content-Sha256' : hashed_payload
            }
        # For POST request, no canonical query string
        canonical_query_string = ''
        # Construct full query string dictionary
        gmtime = time.gmtime()
        datestamp = time.strftime('%Y%m%d', gmtime)
        headers['X-Amz-Date'] = time.strftime('%Y%m%dT%H%M%SZ', gmtime)
        canonical_headers = '\n'.join([(key.lower().strip() + ':' 
                                        + value.strip())
                                        for key, value 
                                        in sorted(headers.items())]) \
                                     + '\n'
        signed_headers \
            = ';'.join([header.lower() for header in sorted(headers)]) 
        canonical_request = \
            '\n'.join(['POST', self.canonical_uri, canonical_query_string,
                        canonical_headers, signed_headers, hashed_payload])
        hashed_canonical_request = hashlib.sha256(
                                            canonical_request
                                        ).hexdigest()
        credential_scope = '/'.join([datestamp, self.region, self.service,
                                        'aws4_request'])
        string_to_sign = '\n'.join([self.algorithm, headers['X-Amz-Date'],
                                    credential_scope,
                                    hashed_canonical_request])
        signature = aws_signature(
                            string_to_sign,
                            self._aws_secret_access_key,
                            datestamp=datestamp,
                            region=self.region,
                            service=self.service
                        )
        headers['Authorization'] = ('%s Credential=%s/%s/%s/%s/aws4_request, '
                'SignedHeaders=%s, Signature=%s'
            ) % (self.algorithm, self._aws_access_key, datestamp, self.region,
                    self.service, signed_headers, signature)
        request = urllib2.Request(''.join([self.base_prefix, self.host_header, 
                                        self.canonical_uri]), 
                                    headers=headers,
                                    data=payload)
        return urlopen_with_retry(request)

def run_job_flow(branding, json_config, force, region='us-east-1', 
                    no_browser=False, s3cmd_exe='s3cmd',
                    s3cmd_config_file=None):
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
        s3cmd_exe: path to s3cmd executable
        s3cmd_config_file: s3cmd configuration file. File should have lines
            with access_key and secret_key. If None, default s3cmd
            configuration is used.

        No return value.
    """
    iface = dp_iface.DooplicityInterface(
                        branding=branding,
                        opener='Started job flow submission script on {time}.'
                    )
    try:
        '''Should ultimately check to see which regions are available to user
        w/ API.'''
        if region not in _aws_regions:
            raise RuntimeError('Region (--region) is invalid. Region must be '
                               'one of {%s}.' % 
                               ', '.join(['"%s"' % a_region for a_region in 
                                            _aws_regions]))
        # Use default s3cmd configuration file if none is specified
        if s3cmd_config_file is not None:
            aws_access_key, aws_secret_access_key \
                = keys_from_config_file(s3cmd_config_file)
            s3_ansible = S3Ansible(s3cmd_exe=s3cmd_exe,
                                    config=s3cmd_config_file)
        elif os.name == 'nt' and os.getenv('USERPROFILE'):
            s3cmd_config_file = os.path.join(os.getenv('USERPROFILE'),
                    'Application Data', 's3cmd.ini'
                )
            aws_access_key, aws_secret_access_key \
                = keys_from_config_file(s3cmd_config_file)
            s3_ansible = S3Ansible(s3cmd_exe=s3cmd_exe,
                                    config=s3cmd_config_file)
        elif os.getenv('HOME'):
            s3cmd_config_file = os.path.join(os.getenv('HOME'), '.s3cfg')
            aws_access_key, aws_secret_access_key \
                = keys_from_config_file(s3cmd_config_file)
            s3_ansible = S3Ansible()
        else:
            raise RuntimeError('Dooplicity failed to find a default s3cmd '
                               'configuration file, and no configuration file '
                               'was specified in its place. If s3cmd is '
                               'indeed installed properly, look for the file '
                               '"s3cmd.ini" on your filesystem in Windows '
                               'or the file ".s3cfg" in Mac OS/Unix. Specify '
                               'the full path of the configuration file '
                               'when running this program again using the '
                               '--s3cfg option.')
        aws_ansible = AWSAnsible(aws_access_key, aws_secret_access_key,
                                    region=region)
        # (In)security?
        del aws_access_key
        del aws_secret_access_key
        # Serialize relevant parts of JSON configuration
        with open(json_config) as json_stream:
            full_payload = json.load(json_stream)
        try:
            job_flow = full_payload['Steps']
        except Exception:
            raise RuntimeError(
                    'JSON file not in proper format. Ensure that the JSON '
                    'object has a Steps key.'
                )
        step_count = len(job_flow)
        steps = OrderedDict()
        # Check steps for requred data
        required_data = set(['input', 'output', 'mapper', 'reducer'])
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
        except KeyError, IndexError:
            iface.fail(
                    'JSON file not in proper format. Ensure '
                    'that each StepConfig object has a HadoopJarStep '
                    'object with an Args array and a Name string.'
                )
            raise
        iface.step('Read job flow from input JSON.')
        '''Check steps for required Hadoop Streaming command-line parameters
        and for whether outputs are writable.'''
        missing_data = defaultdict(list)
        bad_output_data = []
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
                    and dp.Url(step_data['output']).is_s3 \
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
            errors.extend(['Output directory "%s" of step "%s" already '
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
            if dp.Url(step_data['output']).is_s3:
                s3_ansible.remove_dir(steps[step]['output'])
            buckets.add(bucket_from_url(steps[step]['output']))
        for bucket in buckets:
            s3_ansible.create_bucket(bucket)
        iface.step('Set up output directories on S3.')
        job_flow_response = aws_ansible.post_request(json_config)
        json_response = json.load(job_flow_response)
        if 'JobFlowId' not in json_response:
            raise RuntimeError('Job submission failed. Server returned the '
                               'following response: %s.'
                                % json.dumps(json_response, sort_keys=True,
                                    indent=4, separators=(',', ': ')))
        else:
            job_flow_id = json_response['JobFlowId']
            iface.step(('Submitted job flow.\n'
                        '*****Job flow ID is %s .*****\n' % job_flow_id) + 
                        ('*****Submission can be monitored at %s .*****' %
                            _emr_url.format(job_flow_id=job_flow_id)))
            if not no_browser:
                webbrowser.open(_emr_url.format(job_flow_id=job_flow_id),
                                    new=2) # Open in new tab
                iface.step('Opening URL in default browser.')
        iface.done(closer=('Finished job flow submission script '
                           'on {date}. Run time was {length} seconds.')
                )
    except Exception, GeneratorExit:
         # GeneratorExit added just in case this happens on modifying code
        iface.fail()
        raise
    except KeyboardInterrupt, SystemExit:
        if job_flow_id not in locals():
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
    parser.add_argument(
            '-c', '--s3cfg', type=str, required=False, default=None,
            help='s3cmd configuration file. If none is specified, the ' \
                 'default configuration file is used.'
        )
    parser.add_argument(
            '-x', '--s3cmd-exe', type=str, required=False, default='s3cmd',
            help='Path to s3cmd executable. If none is specified, s3cmd ' \
                 'must be in PATH.'
        )
    parser.add_argument(
            '-r', '--region', type=str, required=False, default='us-east-1',
            help='Amazon data center in which to run job.'
        )
    parser.add_argument('-n', '--no-browser', action='store_const', const=True,
            default=False,
            help='Do not automatically open browser after job flow '
                 'submission.')
    dp_iface.add_args(parser)
    args = parser.parse_args(sys.argv[1:])
    run_job_flow(args.branding, args.json_config, args.force, args.region,
                    args.no_browser, args.s3cmd_exe, args.s3cfg)