"""
filemover.py
Part of Dooplicity framework

Utilities for moving files around among local, hdfs, s3, http and ftp
filesystems.
"""

import os
import sys
import subprocess
import time
import random
import threading

def add_args(parser):
    """ Sets up arguments related to moving files around. """
    parser.add_argument(\
        '--s3cfg', metavar='STR', type=str, required=False,
        help='s3cmd configuration file to use (only relevant if some ' \
             'output is being pushed to S3)')
    parser.add_argument(\
        '--acl-public', action='store_const', const=True, default=False,
        help='Make files uploaded to S3 publicly-readable (only relevant if '
             'some output is being pushed to S3)')

class CommandThread(threading.Thread):
    """ Runs a command on a separate thread. """
    def __init__(self, command_list):
        super(CommandThread, self).__init__()
        self.command_list = command_list
        self.process_return = None
    def run(self):
        self.process_return \
            = subprocess.Popen(self.command_list, stdout=sys.stderr).wait()

class FileMover:
    """ Responsible for details on how to move files to and from URLs. """
    
    def __init__(self, args=None, s3cmd_exe='s3cmd', s3cred=None,
                    s3public=False):
        if args is not None:
            self.s3cred, self.s3public = args.s3cfg, args.acl_public
        else:
            self.s3cred, self.s3public = s3cred, s3public
        self.s3cmd_exe = s3cmd_exe
    
    def put(self, filename, url):
        """ Uploads a local file to a URL.

            filename: path to file to upload
            url: URL to which file should be uploaded. Can be directory name
                + '/' or actual file.

            No return value.
        """
        assert os.path.exists(filename)
        if url.is_s3:
            command_list = [self.s3cmd_exe]
            if self.s3cred is not None:
                command_list.append('-c')
                command_list.append(self.s3cred)
            command_list.append('sync')
            if self.s3public:
                command_list.append("--acl-public")
            command_list.append(filename)
            command_list.append(url.to_nonnative_url())
        elif url.is_curlable:
            raise RuntimeError('Can\'t upload to http/ftp URLs.')
        elif url.is_local:
            try:
                os.makedirs(url.to_url())
            except OSError:
                # Directory exists
                pass
            command_list = ['cp', filename, url.to_url()]
        else:
            command_list = ['hadoop', 'fs', '-put']
            command_list.append(filename)
            command_list.append('/'.join([url.to_url(),
                                            os.path.basename(filename)]))
        command = ' '.join(command_list)
        exit_level = subprocess.Popen(command_list, stdout=sys.stderr).wait()
        if exit_level > 0:
            raise RuntimeError('Non-zero exitlevel %d from push command "%s".'
                               % (exit_level, command))

    def exists(self, url):
        if url.is_local:
            return os.path.exists(url.to_url())
        elif url.is_curlable:
            curl_process = subprocess.Popen(['curl', '--head', url.to_url()],
                                shell=True, bufsize=-1, stderr=subprocess.PIPE,
                                stdout=subprocess.PIPE)
            curl_err = curl_process.stderr.read()
            curl_process.wait()
            if 'resolve host' in curl_err:
                return False
            return True
        elif url.is_s3:
            command_list = [self.s3cmd_exe]
            if self.s3cred is not None:
                command_list.append('-c')
                command_list.append(self.s3cred)
            command_list.extend(['ls', url])
            s3cmd_process = subprocess.Popen(command_list, 
                                                stdout=subprocess.PIPE,
                                                stderr=subprocess.PIPE)
            url_string = url.to_url()
            for line in s3cmd_process.stdout:
                if (url_string + '\n') in line or (url_string + '/\n') in line:
                    return True
            s3cmd_process.wait()
            return False
        else:
            return False
    
    def get(self, url, dest="."):
        """ Get a file to local directory.

            url: Remote location of file.
            dest: Local destination of file.

            No return value.
        """
        if url.is_s3:
            command_list = ['s3cmd']
            if self.s3cred is not None:
                command_list.append("-c")
                command_list.append(self.s3cred)
            command_list.append("get")
            command_list.append(url.to_nonnative_url())
            command_list.append(dest)
            command = ' '.join(command_list)
            exit_level \
                = subprocess.Popen(command_list, stdout=sys.stderr).wait()
            if exit_level > 0:
                raise RuntimeError('Non-zero exitlevel %d from s3cmd '
                                   'get command "%s"' % (exit_level, command))
        elif url.is_curlable:
            oldp = os.getcwd()
            os.chdir(dest)
            command_list = ['curl', '-O', '--connect-timeout', '60']
            command_list.append(url.to_url())
            command = ' '.join(command_list)
            while True:
                curl_thread = CommandThread(command_list)
                curl_thread.start()
                while curl_thread.is_alive():
                    print >>sys.stderr, 'reporter:status:alive'
                    sys.stderr.flush()
                    time.sleep(60)
                if curl_thread.process_return > 89:
                    '''If the exit code is greater than the highest-documented
                    curl exit code, there was a timeout.'''
                    print >>sys.stderr, 'Too many simultaneous connections; ' \
                                        'restarting in 10 s.'
                    time.sleep(10)
                else:
                    break
            os.chdir(oldp)
            if curl_thread.process_return > 0:
                raise RuntimeError(('Nonzero exitlevel %d from curl command '
                                    '"%s"') 
                                        % (curl_thread.exit_level, command))
        elif url.is_local:
            command_list = ['cp', url.to_url(), dest]
        else:
            command_list = ["hadoop", "fs", "-get"]
            command_list.append(url.to_url())
            command_list.append(dest)
            command = ' '.join(command_list)
            exit_level \
                = subprocess.Popen(command_list, stdout=sys.stderr).wait()
            if exit_level > 0:
                raise RuntimeError('Nonzero exitlevel %d from hadoop fs '
                                   '-get command "%s"' % (exit_level, command))
