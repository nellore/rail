#!/usr/bin/env python
"""
dooplicity (working name)

Helpful tools for working with Hadoop in Python; may evolve into barebones
framework.

Copyright (c) 2014 Abhi Nellore and Ben Langmead. License TBD.
"""

class xstream:
    """ Permits Pythonic iteration through Hadoop input streams.

        All iterators are implemented as generators.

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
        hadoop_stream: where to find input lines
        key_fields: the first "key_fields" fields from an input line are
            considered the key denoting a partition
        separator: delimiter separating fields from each input line
        skip_duplicates: skip any duplicate lines that may follow a line
    """
    def __init__(self, 
                    hadoop_stream,
                    key_fields=1,
                    separator='\t',
                    skip_duplicates=False):
        if not isinstance(key_fields, int) or key_fields < 1:
            raise RuntimeError('key_fields must be an integer >= 1; value of ' 
                                + str(key_fields) + 'was entered')
        if not isinstance(separator, str) or not len(separator):
            raise RuntimeError('separator must be a string of length '
                               'at least 1')
        self._hadoop_stream = hadoop_stream
        self._key_fields = key_fields
        self._separator = separator
        self._skip_duplicates = skip_duplicates
        self._line = hadoop_stream.readline()
        if self._line:
            self._tokens = self._line.strip().split(separator)
            self._key = tuple(self._tokens[:key_fields])
            self._value = tuple(self._tokens[key_fields:])

    @property
    def key(self):
        try:
            return self._key
        except NameError:
            return None

    @property
    def value(self):
        try:
            return self._value
        except NameError:
            return None

    def _partition(self):
        yield self._value
        self._last_line = self._line
        self._line = self._hadoop_stream.readline()
        if self._skip_duplicates:
            while self._last_line == self._line:
                self._last_line = self._line
                self._line = self._hadoop_stream.readline()
        self._tokens = self._line.strip().split(self._separator)
        self._next_key = tuple(self._tokens[:self._key_fields])
        while self._next_key == self._key:
            self._value = tuple(self._tokens[self._key_fields:])
            yield tuple(self._tokens[self._key_fields:])
            self._last_line = self._line
            self._line = self._hadoop_stream.readline()
            if self._skip_duplicates:
                while self._last_line == self._line:
                    self._last_line = self._line
                    self._line = self._hadoop_stream.readline()
            self._tokens = self._line.strip().split(self._separator)
            self._next_key = tuple(self._tokens[:self._key_fields])
            self._value = tuple(self._tokens[self._key_fields:])
        self._key = self._next_key
        return

    def __iter__(self):
        while self._line:
            self._last_line = None
            yield self._key, self._partition()
            if self._last_line is None:
                raise RuntimeError('Cannot yield next partition without '
                                   'first iterating through values in current '
                                   'partition.')
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

        def test_partitioning(self):
            """ Fails if input data isn't partitioned properly. """
            with open(self.input_file, 'w') as input_stream:
                # Create some fake data with two key fields
                input_stream.write(
                        'chr1\t1\ta\t20\t90\n'
                        'chr1\t1\ti\t10\t50\n'
                        'chr1\t1\ti\t30\t70\n'
                        'chr1\t1\ti\t75\t101\n'
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
                        ('chr1', '2') : [
                                            ('i', '90', '1300'),
                                            ('i', '91', '101')
                                        ]
                    }
                )

        def test_runtime_errors(self):
            """ Fails if RuntimeError isn't raised. """
            with open(self.input_file, 'w') as input_stream:
                # Create some fake data with two key fields
                input_stream.write(
                        'chr1\t1\ta\t20\t90\n'
                        'chr1\t1\ti\t10\t50\n'
                        'chr1\t1\ti\t30\t70\n'
                        'chr1\t1\ti\t75\t101\n'
                        'chr1\t2\ti\t90\t1300\n'
                        'chr1\t2\ti\t91\t101\n'
                    )
            with open(self.input_file) as input_stream:
                partitions = {}
                with self.assertRaises(RuntimeError):
                    for _, _ in xstream(input_stream, 2):
                        pass
            with open(self.input_file) as input_stream:
                with self.assertRaises(RuntimeError):
                    myxstream = xstream(input_stream, '2')
            with open(self.input_file) as input_stream:
                with self.assertRaises(RuntimeError):
                    myxstream = xstream(input_stream, 2, '')

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