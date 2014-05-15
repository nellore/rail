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
import base64
import hashlib
import hmac
import urllib
import urllib2
import json

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
    if url[0] == '/':
        start_index = 1
    elif url[:6] == 's3n://':
        start_index = 6
    elif url[:5] == 's3://':
        start_index = 5
    while url[start_index] == '/':
        start_index += 1
    return url[start_index:url[start_index].index('/')]

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
                            stderr=subprocess.PIPE
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
                                stderr=subprocess.PIPE
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

            path: directory name on #3, including bucket

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
                                stdout=subprocess.PIPE
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

class AWSAnsible:
    """ Interacts with AWS via GET requests. See AWS API documentation. 

        Only Signature Version 2 is supported right now; it works with EMR
        but not all services.

        http://stackoverflow.com/questions/1088715/
        how-to-sign-amazon-web-service-requests-from-the-python-app-engine

        helped make this class.
    """
    def __init__(self, aws_access_key, aws_secret_access_key,
                    base_uri='https://elasticmapreduce.amazonaws.com'):
        self._aws_secret_access_key = aws_secret_access_key
        self._most_authentication_parameters = {
            'AWSAccessKeyId' : aws_access_key,
            'SignatureVersion' : '2',
            'SignatureMethod' : 'HmacSHA256'
        }
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
        self.host_header = base_suffix[0]
        try:
            self.path_segments = '/' + '/'.join([urllib.quote_plus(segment) 
                for segment in base_suffix[1:]])
        except IndexError:
            # No segments
            self.path_segments = '/'

    def get_request(json_object):
        """ Submits JSON payload to AWS using Signature Version 2.

            json_object: JSON object; structured as a tree of dictionaries,
                lists, and strings

            Return value: file-like object of type urllib2.urlopen with
                response
        """
        # Construct full query string dictionary
        query_dict = dict(self._authentication_params.items()
                            + aws_params_from_json(json_object).items())
        query_dict['Timestamp'] \
            = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
        # Sort by key
        keys = query_dict.keys()
        keys.sort()
        query_string = urllib.urlencode(zip(keys, map(query_dict.get, keys)))
        string_to_sign = "GET\n%s\n%s\n%s" % (self.host_header,
                                                self.path_segments,
                                                query_string)
        # Sign, seal, and deliver
        signature_string = '&Signature=%s' %
            urllib.quote_plus(
                    base64.encodestring(
                        hmac.new(
                                    key=self.,
                                    msg=string_to_sign,
                                    digestmod=hashlib.sha256
                            ).digest()
                    ).strip()
                )
        return urllib2.urlopen(''.join([self.base_prefix, self.host_header, 
                                        self.path_segments, '?',
                                        query_string, signature_string]))