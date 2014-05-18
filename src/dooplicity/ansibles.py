"""
ansibles.py
Part of Dooplicity framework

Tools for interacting with local filesystems, the web, and cloud services.

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

import os
import time
import urllib
import urllib2
# For signing requests as per AWS API
import hashlib
import hmac
import binascii
import math
import subprocess
from version import version as _version
from functools import wraps

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
    """ Permits simple interactions with S3 via the AWS CLI. 

        Another option is to use Boto directly, but it comes with a lot of
        unused classes --- tools never touched by Dooplicity/Rail. Further, my
        preference is to integrate tools with good command-line APIs so
        languages can be mixed and matched. That said, the AWS CLI does 
        depend on Boto.... -AN
    """

    def __init__(self, aws_exe='aws', iface=None):
        """ Initializes S3Ansible.

            aws_exe: path to aws executable. Should be unnecessary on all
                platforms if the AWS CLI is installed properly, but users
                may want to specify a different executable if multiple versions
                are installed on their systems.
            iface: object of type DooplicityInterface to which to write status
                updates or None if no such updates should be written.
        """
        self.aws = aws_exe
        self._osdevnull = open(os.devnull, 'w')
        self.iface = iface

    def create_bucket(self, bucket):
        """ Creates a new bucket only if it doesn't already exist.

            Raises subprocess.CalledProcessError if a bucket already exists.

            bucket: string with bucket name to create; can have an s3:// at
                the beginning or not.

            No return value.
        """
        cleaned_url = clean_url(bucket)
        subprocess.check_call(' '.join([self.aws, 's3 mb', cleaned_url]),
                                bufsize=-1,
                                shell=True,
                                stderr=self._osdevnull,
                                stdout=self._osdevnull)

    def is_dir(self, path):
        """ Checks whether a directory on S3 exists.

            There are no directories on S3, really, so is_dir() actually
            returns True if there are object keys on S3 that begin with the
            directory name plus a slash.

            path: directory name on S3, including bucket

            Return value: True iff directory exists; else False
        """
        cleaned_url = clean_url(path) + '/'
        while cleaned_url[-2] == '/':
            cleaned_url = cleaned_url[:-1]
        # Now the URL has just one trailing slash
        try:
            aws_command = ' '.join([self.aws, 's3 ls', cleaned_url])
            list_out = \
                subprocess.check_output(aws_command,
                                         stderr=self._osdevnull,
                                         shell=True,
                                         bufsize=-1).strip()
        except subprocess.CalledProcessError as e:
            # Bucket doesn't exist
            return False
        if list_out:
            # Even if the directory is empty, aws s3 ls outputs something
            return True
        return False

    def remove_dir(self, path):
        """ Removes a directory on S3.

            Technically, this function deletes all object keys that begin with
            the specified identifier; so it's VERY important to append a 
            slash to the end of a path so other directories aren't deleted.

            path: directory name on S3, including bucket

            No return value.
        """
        cleaned_url = clean_url(path) + '/'
        while cleaned_url[-2] == '/':
            cleaned_url = cleaned_url[:-1]
        # Now the URL has just one trailing slash
        aws_command = ' '.join([self.aws, 's3 rm --recursive', cleaned_url])
        if self.iface:
            rm_process = subprocess.Popen(aws_command,
                                            bufsize=-1,
                                            stdout=subprocess.PIPE,
                                            stderr=self._osdevnull,
                                            shell=True)
            for line in rm_process.stdout:
                try:
                    self.iface.status(line.split('\r')[1].strip('\r\n'))
                except IndexError:
                    self.iface.status(line.strip())
            if rm_process.wait():
                '''Bucket does not exist; OR there's another error, but fail
                silently.'''
                pass
        else:
            try:
                subprocess.check_call(aws_command_list, bufsize=-1,
                                            stdout=self._osdevnull,
                                            stderr=self._osdevnull)
            except subprocess.CalledProcessError:
                # Bucket does not exist
                pass

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
    """ Interacts with AWS via POST requests (so far). See AWS API
        documentation. GET requests are cheaper, and if ever a GET request
        need be implemented, the function aws_params_from_json() in this file
        converts JSON into parameters from a query string.

        Uses Signature Version 4, the latest version as of 5/15/14.
        Important: relies on proper installation of AWS CLI to obtain
        keys.

        http://stackoverflow.com/questions/1088715/
        how-to-sign-amazon-web-service-requests-from-the-python-app-engine

        helped make this class.
    """
    def __init__(self, profile='default', region=None, valid_regions=None,
                    base_uri='https://elasticmapreduce.amazonaws.com'):
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
        try:
            self.canonical_uri = '/' + '/'.join([urllib.quote_plus(segment) 
                for segment in split_base_suffix[1:]])
        except IndexError:
            # No segments
            self.canonical_uri = '/'
        self.algorithm = 'AWS4-HMAC-SHA256'
        self._aws_access_key_id = None
        self._aws_secret_access_key = None
        self.region = None
        if profile == 'default':
            # Search environment variables for keys first if profile is default
            try:
                self._aws_access_key_id = os.environ['AWS_ACCESS_KEY_ID']
                self._aws_secret_access_key \
                    = os.environ['AWS_SECRET_ACCESS_KEY']
                to_search = None
            except KeyError:
                to_search = '[default]'
            try:
                # Also grab region
                self.region = os.environ['AWS_DEFAULT_REGION']
            except KeyError:
                pass
        else:
            to_search = '[profile ' + profile + ']'
        # Now search AWS CLI config file for the right profile
        if to_search is not None:
            config_file = os.path.join(os.environ['HOME'], '.aws', 'config')
            try:
                with open(config_file) as config_stream:
                    for line in config_stream:
                        if line.strip() == to_search:
                            break
                    for line in config_stream:
                        tokens = [token.strip() for token in line.split('=')]
                        if tokens[0] == 'region':
                            self.region = tokens[1]
                        elif tokens[0] == 'aws_access_key_id':
                            self._aws_access_key_id = tokens[1]
                        elif tokens[0] == 'aws_secret_access_key':
                            self._aws_secret_access_key = tokens[1]
                        else:
                            line = line.strip()
                            if line[0] == '[' and line[-1] == ']':
                                # Break on start of new profile
                                break
            except IOError:
                raise RuntimeError(('No valid AWS CLI configuration found. '
                                    'Make sure the AWS CLI is installed '
                                    'properly and that one of the following '
                                    'is true:\n1) The environment variables '
                                    '"AWS_ACCESS_KEY_ID" and '
                                    '"AWS_SECRET_ACCESS_KEY" are set to '
                                    'the desired AWS access key ID and '
                                    'secret access key, respectively, or\n'
                                    '2) A file ".aws/config" exists in your '
                                    'home directory with the right settings. '
                                    'To set this file up, run "aws --config" '
                                    'after installing the CLI.'))
        if region is not None:
            # Always override region
            self.region = region
        if self.region is None:
            # If after all this, the region is None, default to us-east-1
            self.region = 'us-east-1'
        # Check for region validity
        if valid_regions is not None:
            if self.region not in valid_regions:
                raise RuntimeError(('Region "%s" is invalid for service '
                                    '"%s". Region must be one of {%s}.') % 
                                        (self.region, self.service,
                                        ', '.join(['"%s"' % a_region
                                                   for a_region
                                                   in valid_regions]))
                                    )
        self.host_header = '.'.join([self.region, split_base_suffix[0]])

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
            ) % (self.algorithm, self._aws_access_key_id, datestamp, self.region,
                    self.service, signed_headers, signature)
        request = urllib2.Request(''.join([self.base_prefix, self.host_header, 
                                        self.canonical_uri]), 
                                    headers=headers,
                                    data=payload)
        return urlopen_with_retry(request)

class Url:
    def __init__(self, url):
        """ Uses prefix to determine type of URL.

            Prefix is part of URL before colon, if it's present.

            url: URL string
        """
        if ':' in url:
            colon_index = url.index(':')
            prefix = url[:colon_index].lower()
            if prefix[:3] == 's3n':
                self.type = 's3n'
            elif prefix[:2] == 's3':
                self.type = "s3"
            elif prefix[:4] == 'hdfs':
                self.type = 'hdfs'
            elif prefix[:4] == 'http':
                self.type = 'http'
            elif prefix[:4] == 'ftp':
                self.type = 'ftp'
            elif prefix[:5] == 'local':
                self.type = 'local'
            else:
                raise RuntimeError(('Unrecognized URL %s; it\'s not S3, HDFS, '
                                   'HTTP, FTP, or local.') % url)
            self.suffix = url[colon_index+1:]
        else:
            self.type = 'local'
            self.suffix = url
        self.is_s3 = self.type[:2] == 's3'
        self.is_curlable = self.type in ['ftp', 'http', 'https']
        self.is_local = self.type == 'local'

    def to_url(self):
        """ Returns URL string: an absolute path if local or an URL.

            Return value: URL string
        """
        if self.type == 'local':
            absolute_path = os.path.abspath(self.suffix)
            return (absolute_path + '/') if self.suffix[-1] == '/' \
                else absolute_path
        else:
            return self.type + ':' + self.suffix

    def plus(self, file_or_subdirectory):
        """ Returns a new URL + file_or_subdirectory.

            file_or_subdirectory: file or subdirectory of URL

            Return value: Url object with file_or_directory tacked on
        """
        original_url = self.to_url()
        return Url(os.path.join(original_url, file_or_subdirectory))

    def upper_url(self):
        """ Uppercases an URL's prefix or gives a local URL's absolute path.

            This is Useful for hiding protocol names from Elastic MapReduce so
            it doesn't mistake a URL passed as a mapper argument as an input
            URL.

            Return value: transformed URL string
        """
        if self.type == 'local':
            absolute_path = os.path.abspath(self.suffix)
            return (absolute_path + '/') if self.suffix[-1] == '/' \
                else absolute_path
        else:
            return self.type.upper() + ':' + self.suffix

    def to_nonnative_url(self):
        """ Converts s3n:// URLs to s3:// URLs 

            Return value: converted URL
        """
        if self.type[:2] == 's3':
            return 's3:' + self.rest
        else:
            return self.type + ':' + self.rest

    def to_native_url(self):
        """ Converts s3:// URLs to s3n:// URLs 

            Return value: converted URL
        """
        if self.type[:2] == 's3':
            return 's3n:' + self.suffix
        else:
            return self.type + ':' + self.suffix
