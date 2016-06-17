#!/usr/bin/env python

"""
counters.py
Part of Dooplicity framework

Abstractions for maintaining and printing counter information.

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

import sys
from collections import defaultdict


class Counter(object):
    """
    Keeps track of counters and flushes them.
    """

    def __init__(self, group, output_fh=sys.stderr, report_style='hadoop'):
        self.counts_flushable = defaultdict(int)
        self.counts = defaultdict(int)
        self.group = group
        self.report_style = report_style
        self.output_fh = output_fh

    def add(self, counter, amt):
        """ Add given amount to given group """
        self.counts[counter] += amt
        self.counts_flushable[counter] += amt

    def get(self, counter):
        """ Return total count associated with group """
        return self.counts[counter]

    def get_since_last_flush(self, counter):
        """ Return count since last call to flush associated with group """
        return self.counts_flushable[counter]

    def flush(self):
        """
        Probably called by atexit or in a keep-alive loop.
        Hadoop streaming counter format: reporter:counter:<group>,<counter>,<amount>
        """
        for k, v in sorted(self.counts_flushable.items()):
            if v == 0:
                continue
            if self.report_style == 'hadoop':
                self.output_fh.write('reporter:counter:')
            self.output_fh.write('%s,%s,%d\n' % (self.group, k, v))
        self.counts_flushable = defaultdict(int)


if __name__ == '__main__':
    import unittest
    from StringIO import StringIO

    if '--test' in sys.argv:

        class TestCounter(unittest.TestCase):

            def test_small_1(self):
                c = Counter('test_sm_1', StringIO())
                c.add('blah', 2)
                c.add('blah', 0)
                self.assertEqual(c.get('blah'), 2)
                self.assertEqual(c.get_since_last_flush('blah'), 2)
                self.assertEqual(c.get('unknown1'), 0)
                self.assertEqual(c.get_since_last_flush('unknown2'), 0)
                c.flush()
                c.add('blah', -20)
                self.assertEqual(c.get('blah'), -18)
                self.assertEqual(c.get_since_last_flush('blah'), -20)

            def test_small_2(self):
                outp = StringIO()
                c = Counter('test_sm_2', outp)
                c.add('ABC', 1)
                c.add('XYZ', 1)
                c.add('XYZ', 1)
                c.add('XYZ', 1)
                c.flush()
                lines = outp.getvalue().rstrip().split('\n')
                self.assertEqual(2, len(lines))
                self.assertEqual(lines[0], "reporter:counter:test_sm_2,ABC,1")
                self.assertEqual(lines[1], "reporter:counter:test_sm_2,XYZ,3")

            def test_small_3(self):
                outp = StringIO()
                c = Counter('test_sm_3', outp)
                c.flush()
                self.assertEqual(0, len(outp.getvalue()))

        unittest.main(argv=[sys.argv[0]])
        sys.exit(0)
