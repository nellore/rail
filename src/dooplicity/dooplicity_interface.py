#!/usr/bin/env python
"""
dooplicity_interface.py

Could later be updated to use curses; however, curses can be messy and give
rise to compatibility issues.
"""
import threading
from dooplicity_version import version
import sys
import time
import tempfile
import os
import json

class style:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'

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
    
    def run(self):
        while True:
            if self.stop.isSet():
                break
            m, s = divmod(time.time() - self._start_time, 60)
            h, m = divmod(m, 60)
            # Edit just the time line of the header
            sys.stdout.write('\r\x1b[K%02dh:%02dm:%02ds | %s' 
                                % (h, m, s, self.message))
            sys.stdout.flush()
            time.sleep(0.2)

class DooplicityInterface:
    def __init__(self, branding=None):
        sys.stdout.write('\n')
        sys.stdout.flush()
        try:
            with open(branding) as branding_stream:
                print branding_stream.readlines()
            print 'Powered by Dooplicity v%s' % version
        except TypeError:
            # No branding
            print 'Dooplicity v%s' % version
        self._start_time = time.time()
        self._date_format = '%A, %b %d, %Y at %I:%M:%S %p'
        if sys.stderr.isatty():
            self._write_streams = [sys.stderr]
        else:
            self._write_streams = [sys.stderr, sys.stdout]
        for output_stream in self._write_streams:
            print >>output_stream, 'Started job on %s.' % time.strftime(
                                    self._date_format,
                                    time.localtime(self._start_time)
                                )
        print >>sys.stdout, style.RED + '\n~.oOo.~\n'
        sys.stdout.flush()
        self._update_thread = UpdateThread(self._start_time)

    def step(self, message):
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
             opener='**Dooplicity encountered errors.**'):
        """ Optionally prints a fail message and the time a job took.

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
        print >>sys.stderr, opener
        if message:
            print >>sys.stderr, message
        end_time = time.time()
        for output_stream in self._write_streams:
            print >>output_stream, 'Job FAILED on %s. Run time was %.03f ' \
                'seconds.' % (time.strftime(self._date_format,
                              time.localtime(end_time)),
                              end_time - self._start_time)
            output_stream.flush()
        if steps:
            temp_dir = tempfile.mkdtemp()
            temp_json_file = os.path.join(temp_dir, 'temp.json')
            with open(temp_json_file, 'w') as json_stream:
                json.dump(steps, json_stream)
                for output_stream in self._write_streams:
                    print >>output_stream, 'To start this job from where it ' \
                                           'left off after errors are ' \
                                           'resolved, run:'
                    print >>output_stream, '%s -j %s -f' \
                        % (os.path.abspath(sys.argv[0]), temp_json_file)
                    output_stream.flush()

    def done(self, message=''):
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
            print >>sys.stderr, message
        end_time = time.time()
        for output_stream in self._write_streams:
            print >>output_stream, ('Job FINISHED on %s. Run time was '
                                    '%.03f seconds.') \
                                    % (time.strftime(self._date_format,
                                      time.localtime(end_time)),
                                      end_time - self._start_time)
            output_stream.flush()