#!/usr/bin/env python
"""
dooplicity_interface.py
Part of Dooplicity framework

Contains Dooplicity console interface class used by local mode, Hadoop mode,
and EMR mode. Could later be updated to use curses; however, curses can be
messy and give rise to compatibility issues.

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

import threading
from dooplicity_version import version
import sys
import time
import tempfile
import os
import json
import itertools

def inflected(number, word, es=False):
    """ Returns string with word in appropriate form.

        number: an integer
        word: singular form of word
        es: whether the plural form takes an 'es'

        Return value: string with number and word
    """
    if number == 1:
        return str(number) + ' ' + word
    if es:
        return str(number) + ' ' + word + 'es'
    return str(number) + ' ' + word + 's'

class UpdateThread(threading.Thread):
    """ Class for updating status/timer. """
    
    def __init__(self, start_time, message=''):
        super(UpdateThread, self).__init__()
        self._start_time = start_time
        self.stop = threading.Event()
        self.message = message
        self.cycle_chars = ['.', 'o', 'O', '@', '*']
        self.daemon = True
    
    def run(self):
        progress_char_gen = itertools.cycle(self.cycle_chars)
        while True:
            if self.stop.isSet():
                break
            m, s = divmod(time.time() - self._start_time, 60)
            h, m = divmod(m, 60)
            # Edit just the time line of the header
            progress = next(progress_char_gen)
            sys.stdout.write('\r\x1b[K%02dh:%02dm:%02ds %s %s' 
                                % (h, m, s, progress, self.message))
            sys.stdout.flush()
            time.sleep(.04)

class DooplicityInterface:
    """ Encapsulates methods for writing status updates to console. """

    def __init__(self, branding=None):
        # Disable line-wrapping
        sys.stdout.write('\n')
        sys.stdout.flush()
        try:
            with open(branding) as branding_stream:
                sys.stdout.write(branding_stream.read())
        except TypeError:
            # No branding
            print 'Dooplicity v%s' % version
        self._start_time = time.time()
        self._date_format = '%A, %b %d, %Y at %I:%M:%S %p %Z'
        if sys.stderr.isatty():
            self._write_streams = [sys.stderr]
        else:
            self._write_streams = [sys.stderr, sys.stdout]
        for output_stream in self._write_streams:
            print >>output_stream, 'Started job flow on %s.' % time.strftime(
                                    self._date_format,
                                    time.localtime(self._start_time)
                                )
        print '\n~.oOo.~\n'
        sys.stdout.flush()
        self._update_thread = UpdateThread(self._start_time)

    def step(self, message):
        """ Writes a step start/finish message to the console.

            message: string

            No return value.
        """
        # Pause update_thread
        self._update_thread.stop.set()
        try:
            self._update_thread.join()
        except RuntimeError:
            # Thread hasn't started
            pass
        # Clear a line
        sys.stdout.write('\r\x1b[K')
        sys.stdout.flush()
        m, s = divmod(time.time() - self._start_time, 60)
        h, m = divmod(m, 60)
        for output_stream in self._write_streams:
            print >>output_stream, '%02dh:%02dm:%02ds | %s' \
                                    % (h, m, s, message)
            output_stream.flush()
        # Restart self._update_thread
        self._update_thread = UpdateThread(self._start_time)
        self._update_thread.start()

    def status(self, message):
        """ Changes status message next to timer on console.

            message: string

            No return valuee.
        """
        self._update_thread.stop.set()
        try:
            self._update_thread.join()
        except RuntimeError:
            # Thread hasn't started
            pass
        self._update_thread = UpdateThread(self._start_time)
        self._update_thread.message = message
        self._update_thread.start()

    def fail(self, message='', steps=[],
             opener='*****Errors encountered*****'):
        """ Optionally prints a fail message and run time took on failure.

            Also outputs command for resuming job if one or more steps
            completed.

            opener: string containing opening message
            message: string containing error message
            steps: JSON array of residual steps in pipeline

            No return value.
        """
        # Terminate update thread
        self._update_thread.stop.set()
        try:
            self._update_thread.join()
        except RuntimeError:
            # Thread hasn't started
            pass
        # Clear a line
        sys.stdout.write('\r\x1b[K')
        sys.stdout.flush()
        for output_stream in self._write_streams:
            print >>output_stream, opener
        if message:
            for output_stream in self._write_streams:
                print >>output_stream, message
        end_time = time.time()
        for output_stream in self._write_streams:
            print >>output_stream, 'Job flow failed on %s. Run time was ' \
                '%.03f seconds.' % (time.strftime(self._date_format,
                                    time.localtime(end_time)),
                                    end_time - self._start_time)
            output_stream.flush()
        if steps:
            temp_dir = tempfile.mkdtemp()
            temp_json_file = os.path.join(temp_dir, 'temp.json')
            with open(temp_json_file, 'w') as json_stream:
                json.dump(steps, json_stream)
                for output_stream in self._write_streams:
                    print >>output_stream, 'To start this job flow from where '
                                           'it left off, run:'
                    print >>output_stream, '%s %s -j %s -f' \
                        % (sys.executable, os.path.abspath(sys.argv[0]),
                           temp_json_file)
                    output_stream.flush()

    def done(self, message=''):
        """ Writes a message on completion of the job flow.

            message: string

            No return valuee.
        """
        # Terminate update thread
        self._update_thread.stop.set()
        try:
            self._update_thread.join()
        except RuntimeError:
            # Thread wasn't started
            pass
        # Clear a line
        sys.stdout.write('\r\x1b[K')
        sys.stdout.flush()
        if message:
            for output_stream in self._write_streams:
                print >>output_stream, message
        end_time = time.time()
        print '\n~.oOo.~\n'
        sys.stdout.flush()
        for output_stream in self._write_streams:
            print >>output_stream, ('Finished job flow on %s. Run time was '
                                    '%.03f seconds.') \
                                    % (time.strftime(self._date_format,
                                      time.localtime(end_time)),
                                      end_time - self._start_time)
            output_stream.flush()
        sys.stdout.write('\n')
        sys.stdout.flush()