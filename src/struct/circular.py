'''
circular.py

Circular buffers.
'''

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
        self.advanceTo(self._off + 1)
    
    def advanceTo(self, off):
        """ Advance to given offset """
        assert off >= self._off
        for _ in xrange(0, off - self._off):
            self._list[self._cur] = 0 # Clear this elt
            self._off += 1
            self._cur += 1
            if self._cur >= len(self._list):
                self._cur = 0
    
    def get(self, off=None):
        """ Get count at given offset """
        if off is None: off = self._off
        diff = off - self._off
        assert diff < len(self._list) # make sure we didn't look too far ahead
        i = self._cur + diff
        if i >= len(self._list):
            i -= len(self._list)
        return self._list[self._cur]
    
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
        ret = []
        if en - st + 1 > self._maxlen:
            self._maxlen = en - st + 1
            self._ends.resize(self._maxlen)
        storig = self._lastst
        while st > self._lastst:
            # Advance to the position
            ret.append(self._cov)
            assert self._ends.off() == self._lastst
            self._ends.advance()
            self._lastst += 1
            enamt = self._ends.get()
            assert enamt <= self._cov, "%d, %d" % (enamt, self._cov)
            self._cov -= enamt
        self._cov += amt
        self._ends.add(amt, en)
        return (storig, ret)
    
    def finalize(self):
        ret = []
        storig = self._lastst
        while self._en > self._lastst:
            # Advance to the position
            ret.append(self._cov)
            assert self._ends.off() == self._lastst
            self._ends.advance()
            self._lastst += 1
            enamt = self._ends.get()
            assert enamt <= self._cov, "%d, %d" % (enamt, self._cov)
            self._cov -= enamt
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
    
    unittest.main()
