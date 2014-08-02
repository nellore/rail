#!/usr/bin/env python
"""
tools.py
Part of Dooplicity framework

Includes a class for iterating through streams easily and a few other tools.

The functions which() and is_exe() was taken from
http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
and is not covered under the license below.

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

from itertools import groupby
import os
import threading

class KeepAlive(threading.Thread):
    """ Writes Hadoop status messages to avert task termination. """
    def __init__(self, status_stream, period=120):
        """
            status_stream: where to write status messages
            period: number of seconds between successive status messages
        """
        super(KeepAlive, self).__init__()
        self.period = period
        self.status_stream = status_stream
        # Kills thread when script is finished
        self.daemon = True

    def run(self):
        import time
        while True:
            self.status_stream.write('\nreporter:status:alive\n')
            self.status_stream.flush()
            time.sleep(self.period)

def is_exe(fpath):
    """ Tests whether a file is executable.

        fpath: path to file

        Return value: True iff file exists and is executable.
    """
    return os.path.exists(fpath) and os.access(fpath, os.X_OK)

def which(program):
    """ Tests whether an executable is in PATH.

        program: executable to search for

        Return value: path to executable or None if the executable is not
            found.
    """
    def ext_candidates(fpath):
        yield fpath
        for ext in os.environ.get("PATHEXT", "").split(os.pathsep):
            yield fpath + ext

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            for candidate in ext_candidates(exe_file):
                if is_exe(candidate):
                    return candidate

    return None

def path_join(unix, *args):
    """ Performs UNIX-like os.path.joins on Windows systems if necessary.

        unix: True iff UNIX-like path join should be performed; else False

        Return value: joined path
    """
    args_list = []
    if unix:
        for i in xrange(len(args) - 1):
            try:
                if args[i][-1] != '/':
                    args_list.append(args[i] + '/')
            except IndexError:
                # Empty element
                pass
        args_list.append(args[-1])
        return ''.join(args_list)
    else:
        return os.path.join(*args)

def xopen(gzipped, *args):
    """ Passes args on to the appropriate opener, gzip or regular.

        gzipped: True iff gzip.open() should be used to open rather than
            open(); False iff open() should be used; None if input should be
            read and guessed
        *args: unnamed arguments to pass

        Return value: file object
    """
    import gzip
    if gzipped is None:
        with open(args[0], 'rb') as binary_input_stream:
            # Check for magic number
            if binary_input_stream.read(2) == '\x1f\x8b':
                gzipped = True
            else:
                gzipped = False
    if gzipped:
        return gzip.open(*args)
    return open(*args)

class dlist:
    """ List data type that spills to disk if a memlimit is reached.

        Keeping memory usage low can be important in Hadoop, so this class
        is included in Dooplicity.

        Random access is not currently permitted. The list should properly
        be used by appending all elements, then iterating through them to
        read them.
    """
    def __init__(self, limit=1000000):
        """
            limit: maximum number of elements allowed in list before
                spilling to disk
        """
        self.mem_list = []
        self.disk_stream = None
        self.limit = limit

    def __iter__(self):
        """ Iterates through list.

            NOTE THAT seek to beginning of file is performed if some of the
            list is not in memory!
        """
        if self.disk_stream is not None:
            self.disk_stream.flush()
            self.disk_stream.seek(0)
            for line in self.disk_stream:
                yield line.strip()
        else:
            for item in self.mem_list:
                yield item

    def append(self, item):
        """ Appends item to list. Only strings are permitted right now.

            item: string to append
        """
        if type(item) is not str:
            raise TypeError('An item appended to a dlist must be a string.')
        if self.mem_list is not None and len(self.mem_list) < self.limit:
            self.mem_list.append(item)
        elif self.disk_stream is None:
            # Open new temporary file
            import tempfile
            self.disk_stream = tempfile.TemporaryFile()
            print >>self.disk_stream, '\n'.join(self.mem_list)
            self.mem_list = None
            print >>self.disk_stream, item
        else:
            print >>self.disk_stream, item

class xstream:
    """ Permits Pythonic iteration through partitioned/sorted input streams.

        All iterators are implemented as generators. Could have subclassed
        itertools.groupby here; however, implementation of itertools.groupby
        may change from version to version of Python. Implementation is thus
        just based on itertools.groupby from
        https://docs.python.org/2/library/itertools.html .

        Usage: for key, xpartition in xstream(hadoop_stream):
                   for value in xpartition:
                        <code goes here>

        Each of key and value above is a tuple of strings.

        Properties
        -------------
        key: key tuple that denotes current partition; this is an attribute
            of both an xstream. None when no lines have been read yet.
        value: tuple that denotes current value. None when no lines have been
            read yet.

        Init vars
        -------------
        input_stream: where to find input lines
        key_fields: the first "key_fields" fields from an input line are
            considered the key denoting a partition
        separator: delimiter separating fields from each input line
        skip_duplicates: skip any duplicate lines that may follow a line
    """
    @staticmethod
    def stream_iterator(
            input_stream,
            separator='\t',
            skip_duplicates=False
        ):
        if skip_duplicates:
            for line, _ in groupby(input_stream):
                yield tuple(line.strip().split(separator))
        else:
            for line in input_stream:
                yield tuple(line.strip().split(separator))

    def __init__(
            self, 
            input_stream,
            key_fields=1,
            separator='\t',
            skip_duplicates=False
        ):
        self._key_fields = key_fields
        self.it = self.stream_iterator(
                        input_stream,
                        separator=separator,
                        skip_duplicates=skip_duplicates
                    )
        self.tgtkey = self.currkey = self.currvalue = object()

    def __iter__(self):
        return self

    def next(self):
        while self.currkey == self.tgtkey:
            self.currvalue = next(self.it)    # Exit on StopIteration
            self.currkey = self.currvalue[:self._key_fields]
        self.tgtkey = self.currkey
        return self.currkey, self._grouper(self.tgtkey)

    def _grouper(self, tgtkey):
        while self.currkey == tgtkey:
            yield self.currvalue[self._key_fields:]
            self.currvalue = next(self.it)    # Exit on StopIteration
            self.currkey = self.currvalue[:self._key_fields]

if __name__ == '__main__':
    # Run unit tests
    import unittest
    import tempfile
    import os
    import shutil

    class TestXstream(unittest.TestCase):
        """ Tests xstream class. """
        def setUp(self):
            # Set up temporary directory
            self.temp_dir_path = tempfile.mkdtemp()
            self.input_file = os.path.join(self.temp_dir_path, 'hadoop.temp')

        def test_partitioning_1(self):
            """ Fails if input data isn't partitioned properly. """
            with open(self.input_file, 'w') as input_stream:
                # Create some fake data with two key fields
                input_stream.write(
                        'chr1\t1\ta\t20\t90\n'
                        'chr1\t1\ti\t10\t50\n'
                        'chr1\t1\ti\t30\t70\n'
                        'chr1\t1\ti\t75\t101\n'
                        'chr3\t2\ti\t90\t1300\n'
                        'chr1\t2\ti\t90\t1300\n'
                        'chr1\t2\ti\t91\t101\n'
                    )
            with open(self.input_file) as input_stream:
                partitions = {}
                for key, xpartition in xstream(input_stream, 2):
                    partitions[key] = []
                    for value in xpartition:
                        partitions[key].append(value)
            self.assertEqual(partitions, 
                    {
                        ('chr1', '1') : [
                                            ('a', '20', '90'),
                                            ('i', '10', '50'),
                                            ('i', '30', '70'),
                                            ('i', '75', '101')
                                        ],
                        ('chr3', '2') : [
                                            ('i', '90', '1300')
                                        ],
                        ('chr1', '2') : [
                                            ('i', '90', '1300'),
                                            ('i', '91', '101')
                                        ]
                    }
                )

        def test_partitioning_2(self):
            """ Fails if input data isn't partitioned properly. """
            with open(self.input_file, 'w') as input_stream:
                # Create some fake data with two key fields
                input_stream.write(
                        '1\tA\n'
                        '1\tB\n'
                        '1\tC\n'
                        '1\tD\n'
                        '1\tE\n'
                        '1\tF\n'
                        '2\tG\n'
                        '3\tH\n'
                        '3\tI\n'
                        '3\tJ\n'
                    )
            with open(self.input_file) as input_stream:
                partitions = {}
                for key, xpartition in xstream(input_stream, 1):
                    partitions[key] = []
                    for value in xpartition:
                        partitions[key].append(value)
            self.assertEqual(partitions, 
                    {
                        ('1',) : [('A',), ('B',), ('C',),
                                  ('D',), ('E',), ('F',)],
                        ('2',) : [('G',)],
                        ('3',) : [('H',), ('I',), ('J',)]
                    }
                )

        def test_duplicate_line_skipping(self):
            """ Fails if duplicate lines aren't skipped. """
            with open(self.input_file, 'w') as input_stream:
                # Create some fake data with two key fields
                input_stream.write(
                        'chr1\t1\ta\t20\t90\n'
                        'chr1\t1\ta\t20\t90\n'
                        'chr1\t1\ti\t10\t50\n'
                        'chr1\t1\ti\t30\t70\n'
                        'chr1\t1\ti\t30\t70\n'
                        'chr1\t1\ti\t30\t70\n'
                        'chr1\t1\ti\t30\t70\n'
                        'chr1\t1\ti\t75\t101\n'
                        'chr1\t1\ti\t75\t101\n'
                        'chr1\t2\ti\t90\t1300\n'
                        'chr1\t2\ti\t91\t101\n'
                    )
            with open(self.input_file) as input_stream:
                output = [(key, value) for key, xpartition
                            in xstream(input_stream, 2, skip_duplicates=True)
                            for value in xpartition]
            self.assertEqual(output,
                    [(('chr1', '1'), ('a', '20', '90')),
                     (('chr1', '1'), ('i', '10', '50')),
                     (('chr1', '1'), ('i', '30', '70')),
                     (('chr1', '1'), ('i', '75', '101')),
                     (('chr1', '2'), ('i', '90', '1300')),
                     (('chr1', '2'), ('i', '91', '101'))]
                )

        def test_empty_input(self):
            """ Fails if it fails. """
            with open(self.input_file, 'w') as input_stream:
                pass
            with open(self.input_file) as input_stream:
                for key, xpartition in xstream(input_stream, 1):
                    for value in xpartition:
                        pass

        def tearDown(self):
            # Kill temporary directory
            shutil.rmtree(self.temp_dir_path)

    unittest.main()
