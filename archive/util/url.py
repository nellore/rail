"""
url.py
Part of Dooplicity framework

Includes a class for parsing, handling and manipulating URLs that might be on
the local filesystem, in HDFS, or in either of the S3 filesystems
(native or block).

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

class Url(object):
    def __init__(self, url):
        """ Uses prefix to determine type of URL.

            Prefix is part of URL before colon, if it's present.

            url: URL string
        """
        if ':' in url:
            colon_index = s.index(':')
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
        self.is_curlable = self.type in ['ftp', 'http']
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