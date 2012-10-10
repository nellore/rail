'''
circular.py


'''

class CircularCountBuffer(object):
    
    def __init__(self, off, sz):
        self._cur = 0          # cursor
        self._off = off        # offset under the cursor
        self._list = [0] * sz  # elements
    
    def __len__(self):
        return len(self._list)
        
    def resize(self, newsz):
        newlist = [0] * newsz
        cur = self._cur
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
        self.advanceTo(self._off + 1)

    def advanceTo(self, off):
        assert off >= self._off
        for _ in xrange(0, off - self._off):
            self._list[self._cur] = 0 # Clear this elt
            self._off += 1
            self._cur += 1
            if self._cur >= len(self._list):
                self._cur = 0

    def get(self, off=None):
        if off is None: off = self._off
        diff = off - self._off
        assert diff < len(self._list)
        i = self._cur + diff
        if i >= len(self._list):
            i -= len(self._list)
        return self._list[self._cur]

    def inc(self, off=None):
        if off is None: off = self._off
        assert off >= self._off
        diff = off - self._off
        assert diff < len(self._list)
        i = self._cur + diff
        if i >= len(self._list):
            assert off > 0
            i -= len(self._list)
        self._list[i] += 1
    
    def off(self): return self._off

    def cur(self): return self._cur

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

    unittest.main()
