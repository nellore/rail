'''
circular.py

Circular buffers.
'''

import sys

class CircularCountBuffer(object):
    
    """ Buffer for holding integer counts temporarily as we move sequentially
        along a number line.  The buffer has an offset (the offset where we
        begin our movement) and a maximum size (max # positions we'll have to
        keep track of at once). """
    
    def __init__(self, off, sz):
        self._cur = 0          # cursor
        self._off = off        # offset under the cursor
        self._list = [0] * sz  # elements
    
    def __len__(self):
        """ Return size of buffer """
        return len(self._list)
    
    def resize(self, newsz):
        """ Resize buffer """
        newlist = [0] * newsz
        cur = self._cur
        # Copy _list into newlist
        for i in xrange(0, len(self._list)):
            if i >= len(newlist):
                break
            newlist[i] = self._list[cur]
            cur += 1
            if cur >= len(self._list):
                cur = 0
        self._cur = 0
        self._list = newlist
    
    def advance(self):
        """ Advance to next offset """
        return self.advanceTo(self._off + 1)
    
    def advanceTo(self, off):
        """ Advance to given offset """
        assert off >= self._off
        for _ in xrange(0, off - self._off):
            self._list[self._cur] = 0 # Clear this elt
            self._off += 1
            self._cur += 1
            if self._cur >= len(self._list):
                self._cur = 0
        return self._list[self._cur]
    
    def get(self, off=None):
        """ Get count at given offset """
        if off is None: off = self._off
        diff = off - self._off
        assert diff < len(self._list) # make sure we didn't look too far ahead
        i = self._cur + diff
        if i >= len(self._list):
            i -= len(self._list)
        return self._list[self._cur]
    
    def peekMany(self, n):
        ret = [0] * n
        assert self._cur < len(self._list)
        left = min(len(self._list) - self._cur, n)
        if left > 0:
            ret[0:left] = self._list[self._cur:(self._cur+left)]
            assert len(ret) == n
        if left < n:
            # Note, n might be longer than buffer
            ncopy = min(n-left, self._cur)
            assert ncopy < len(self._list)
            ret[left:(left+ncopy)] = self._list[0:ncopy]
            assert len(ret) == n, "left=%d, n=%d, len(ret)=%d" % (left, n, len(ret))
        assert len(ret) == n
        return ret
    
    def inc(self, off=None):
        """ Increment count at given offset """
        self.add(1, off)
    
    def add(self, amt, off=None):
        """ Add to a count at a given offset """
        if off is None: off = self._off
        assert off >= self._off, "Expected given offset (%d) >= circular buffer offset (%d)" % (off, self._off)
        diff = off - self._off
        assert diff < len(self._list) # make sure we didn't look too far ahead
        i = self._cur + diff
        if i >= len(self._list):
            assert off > 0
            i -= len(self._list)
        self._list[i] += amt
    
    def off(self):
        """ Get current offset """
        return self._off
    
    def cur(self):
        """ Get cursor position """
        return self._cur

class CircularCoverageBuffer(object):
    
    """ A buffer that minimizes the amount of memory we have to use to compile
       coverage information over a sorted list of intervals.  Intervals are
       sorted by the starting position along the genome. """
    
    def __init__(self, st, en, maxlen):
        self._st = st
        self._en = en
        self._maxlen = maxlen
        self._ends = CircularCountBuffer(st, maxlen)
        self._lastst = st
        self._cov = 0
    
    def add(self, st, en, amt):
        """ Add an interval; return coverages for fully-resolved positions """
        assert en > st
        assert en > self._st
        assert st <= self._en
        if en - st + 1 > self._maxlen:
            self._maxlen = en - st + 1
            self._ends.resize(self._maxlen)
        storig = self._lastst
        ret = [0] * (st - self._lastst)
        if st > self._lastst:
            self._ends.advance()
            cov = self._ends.peekMany(st - self._lastst)
            for i in xrange(0, st - self._lastst):
                ret[i] = self._cov
                self._cov -= cov[i]
            self._lastst = st
            self._ends.advanceTo(st)
        self._cov += amt
        self._ends.add(amt, en)
        return (storig, ret)
    
    def finalize(self):
        self._ends.advance()
        ret = [0] * (self._en - self._lastst)
        if self._cov > 0:
            cov = self._ends.peekMany(self._en - self._lastst)
            for cur in xrange(0, self._en - self._lastst):
                # Advance to the position
                ret[cur] = self._cov
                self._cov -= cov[cur]
                if self._cov == 0:
                    break
        return (self._lastst, ret)
    
class CircularMultiCoverageBuffer(object):
    
    """ A buffer that minimizes the amount of memory we have to use to compile
       coverage information for multiple samples over a sorted list of
       intervals.  Intervals are sorted by the starting position along the
       genome. """
    
    def __init__(self, nsamps, st, en, maxlen):
        self._st = st
        self._en = en
        self._nsamps = nsamps
        self._maxlen = maxlen
        self._lastst = st
        self._reset()
    
    def _reset(self):
        self._cov = [0] * self._nsamps
        self._ends = [ CircularCountBuffer(self._st, self._maxlen) for _ in xrange(0, self._nsamps) ]
    
    def add(self, sampid, st, en, amt):
        """ Add an interval; return coverages for fully-resolved positions """
        assert sampid < self._nsamps
        assert en > st
        assert en > self._st
        assert st <= self._en
        ret = None
        if en - st + 1 > self._maxlen:
            self._maxlen = en - st + 1
            for e in self._ends:
                e.resize(self._maxlen)
        storig = self._lastst
        while st > self._lastst:
            if ret is None:
                ret = [ [] for _ in xrange(0, self._nsamps) ]
            for i in xrange(0, self._nsamps):
                ret[i].append(self._cov[i])
                assert self._ends[i].off() == self._lastst
                self._ends[i].advance()
                enamt = self._ends[i].get()
                assert enamt <= self._cov[i]
                self._cov[i] -= enamt
            self._lastst += 1
        self._cov[sampid] += amt
        self._ends[sampid].add(amt, en)
        return (storig, ret)
    
    def finalize(self):
        ret = None
        storig = self._lastst
        while self._en > self._lastst:
            if ret is None:
                ret = [ [] for _ in xrange(0, self._nsamps) ]
            for i in xrange(0, self._nsamps):
                ret[i].append(self._cov[i])
                assert self._ends[i].off() == self._lastst
                self._ends[i].advance()
                enamt = self._ends[i].get()
                assert enamt <= self._cov[i]
                self._cov[i] -= enamt
            self._lastst += 1
        return (storig, ret)

if __name__ == '__main__':
    import unittest

    class TestCircularCountBuffer(unittest.TestCase):
        
        def test1(self):
            buf = CircularCountBuffer(27, 10)
            self.assertEqual(10, len(buf))
        
        def test2(self):
            buf = CircularCountBuffer(27, 10)
            buf.inc()
            self.assertEqual(10, len(buf))
            self.assertEqual(1, buf.get())
            buf.inc()
            self.assertEqual(2, buf.get())
        
        def test3(self):
            buf = CircularCountBuffer(27, 10)
            buf.inc()
            self.assertEqual(10, len(buf))
            self.assertEqual(1, buf.get())
            buf.inc()
            self.assertEqual(2, buf.get())
            buf.advance()
            self.assertEqual(0, buf.get())
            buf.inc(29)
            buf.inc(29)
            buf.inc(29)
            self.assertEqual(0, buf.get())
            buf.advance()
            self.assertEqual(3, buf.get())
        
        def test4(self):
            buf = CircularCountBuffer(27, 10)
            self.assertEqual(0, buf.get())
            buf.inc(27)
            self.assertEqual(1, buf.get())
            buf.inc(36)
            buf.advanceTo(35)
            self.assertEqual(0, buf.get())
            buf.advance()
            self.assertEqual(1, buf.get())
            buf.advance()
            self.assertEqual(0, buf.get())
            self.assertEqual(37, buf.off())
            self.assertEqual(0, buf.cur())
    
    class TestCircularCoverageBuffer(unittest.TestCase):
        
        def test1(self):
            buf = CircularCoverageBuffer(0, 100, 10)
            ret1 = buf.add(0, 10, 1)
            ret2 = buf.add(0, 10, 2)
            ret3 = buf.add(4, 10, 2)
            ret4 = buf.add(8, 20, 3)
            self.assertEqual((0, []), ret1)
            self.assertEqual((0, []), ret2)
            self.assertEqual((0, [3, 3, 3, 3]), ret3)
            self.assertEqual((4, [5, 5, 5, 5]), ret4)
        
        def test2(self):
            buf = CircularCoverageBuffer(0, 25, 4)
            ret1 = buf.add(5, 9, 1)
            ret2 = buf.add(6, 7, 2)
            ret3 = buf.add(10, 20, 2)
            ret4 = buf.finalize()
            self.assertEqual((0, [0, 0, 0, 0, 0]), ret1)
            self.assertEqual((5, [1]), ret2)
            self.assertEqual((6, [3, 1, 1, 0]), ret3)
            self.assertEqual((10, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0]), ret4)
        
        def test3(self):
            buf = CircularCoverageBuffer(20, 40, 4)
            ret1 = buf.add(5, 21, 1)
            ret2 = buf.add(6, 21, 2)
            ret3 = buf.add(7, 21, 2) # 5 covering 20
            ret4 = buf.add(10, 25, 2) # 2 covering up to 25
            ret5 = buf.add(30, 50, 10)
            ret6 = buf.finalize()
            self.assertEqual((20, []), ret1)
            self.assertEqual((20, []), ret2)
            self.assertEqual((20, []), ret3)
            self.assertEqual((20, []), ret4)
            self.assertEqual((20, [7, 2, 2, 2, 2, 0, 0, 0, 0, 0]), ret5)
            self.assertEqual((30, [10, 10, 10, 10, 10, 10, 10, 10, 10, 10]), ret6)
        
        def test4(self):
            buf = CircularCoverageBuffer(20, 10000, 100)
            ret1 = buf.add(0, 21, 1)
            ret2 = buf.add(1, 22, 2)
            ret3 = buf.add(2, 23, 2) # 5 covering 20
            ret4 = buf.finalize()
            self.assertEqual((20, []), ret1)
            self.assertEqual((20, []), ret2)
            self.assertEqual((20, []), ret3)
            self.assertEqual(20, ret4[0])
            self.assertEqual([5, 4, 2], ret4[1][0:3])
            self.assertEqual([0] * 9977, ret4[1][3:])
        
        def test5(self):
            import random
            ln = 1000
            max_ival = 250
            buf = CircularCoverageBuffer(0, ln, max_ival)
            ivals = []
            cov = {}
            for _ in xrange(0, 100):
                ival_len = random.randint(40, max_ival)
                off = random.randint(0, ln - ival_len)
                ivals.append((off, off + ival_len, 1))
                for i in xrange(0, ival_len):
                    cov[i+off] = cov.get(i+off, 0) + 1
            ivals.sort()
            cur = 0
            for st, en, amt in ivals:
                bufcov = buf.add(st, en, amt)
                bst, bufcov = bufcov
                self.assertEqual(bst, cur)
                for i in xrange(0, len(bufcov)):
                    covcov = cov.get(cur, 0)
                    self.assertEqual(covcov, bufcov[i])
                    cur += 1
            bufcov = buf.finalize()
            bst, bufcov = bufcov
            covcov = cov.get(cur, 0)
            for i in xrange(0, len(bufcov)):
                covcov = cov.get(cur, 0)
                self.assertEqual(covcov, bufcov[i])
                cur += 1
    
    class TestCircularMultiCoverageBuffer(unittest.TestCase):
        
        def test1(self):
            buf = CircularMultiCoverageBuffer(2, 0, 100, 10)
            ret1 = buf.add(0, 0, 10, 1)
            ret2 = buf.add(0, 0, 10, 2)
            ret3 = buf.add(0, 4, 10, 2)
            ret4 = buf.add(0, 8, 20, 3)
            self.assertEqual((0, None), ret1)
            self.assertEqual((0, None), ret2)
            self.assertEqual((0, [[3, 3, 3, 3], [0, 0, 0, 0]]), ret3)
            self.assertEqual((4, [[5, 5, 5, 5], [0, 0, 0, 0]]), ret4)
        
        def test2(self):
            buf = CircularMultiCoverageBuffer(2, 0, 20, 10)
            ret1 = buf.add(0, 0, 10, 1)
            ret2 = buf.add(0, 0, 10, 2)
            ret3 = buf.add(1, 4, 10, 2)
            ret4 = buf.add(0, 8, 20, 3)
            ret5 = buf.finalize()
            self.assertEqual((0, None), ret1)
            self.assertEqual((0, None), ret2)
            self.assertEqual((0, [[3, 3, 3, 3], [0, 0, 0, 0]]), ret3)
            self.assertEqual((4, [[3, 3, 3, 3], [2, 2, 2, 2]]), ret4)
            self.assertEqual((8, [[6, 6, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3], [2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]), ret5)
        
    unittest.main()
