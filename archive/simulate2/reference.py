""" Encapsulate FASTA files with accompanying index (.fai) """

__author__ = "Ben Langmead"

import re


class ReferenceOOB(Exception):
    """ Out of bounds exception for reference sequences """

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class ReferenceIndexed(object):
    """ Like Reference but uses .fai index files to avoid ever loading
        entire sequences into memory.  Use in Python 'with' block so
        that FASTA filehandles are closed appropriately. """

    __removeWs = re.compile(r'\s+')

    def __init__(self, fafns):
        self.fafhs = {}
        self.faidxs = {}
        self.chr2fh = {}
        self.offset = {}
        self.lens = {}
        self.chars_per_line = {}
        self.bytes_per_line = {}

        for fafn in fafns:
            self.fafhs[fafn] = fh = open(fafn, 'r')
            # Parse the index files
            with open(fafn + '.fai') as idxfh:
                for ln in idxfh:
                    toks = ln.rstrip().split()
                    if len(toks) == 0:
                        continue
                    assert len(toks) == 5
                    ref_id, ln, offset, chars_per_line, bytes_per_line = toks
                    self.chr2fh[ref_id] = fh
                    self.offset[ref_id] = int(offset)  # 0-based
                    self.lens[ref_id] = int(ln)
                    self.chars_per_line[ref_id] = int(chars_per_line)
                    self.bytes_per_line[ref_id] = int(bytes_per_line)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        # Close all the open FASTA files
        for fafh in self.fafhs.itervalues():
            fafh.close()

    def has_name(self, refid):
        return refid in self.offset

    def names(self):
        return self.offset.iterkeys()

    def length(self, refid):
        return self.lens[refid]

    def get(self, refid, start, ln):
        """ Return the specified substring of the reference. """
        if refid not in self.offset:
            raise ReferenceOOB('No such reference sequence name: "%s"' % refid)
        if start + ln > self.lens[refid]:
            raise ReferenceOOB(
                '"%s" has length %d; tried to get [%d, %d)' % (refid, self.lens[refid], start, start + ln))
        fh, offset, chars_per_line, bytes_per_line = self.chr2fh[refid], self.offset[refid], \
                                                     self.chars_per_line[refid], self.bytes_per_line[refid]
        assert bytes_per_line > chars_per_line
        byte_off = offset
        byte_off += (start // chars_per_line) * bytes_per_line
        into = start % chars_per_line
        byte_off += into
        fh.seek(byte_off)
        left = chars_per_line - into
        # Count the number of line breaks interrupting the rest of the
        # string we're trying to read
        if ln < left:
            return fh.read(ln)
        else:
            nbreaks = 1 + (ln - left) // chars_per_line
            res = fh.read(ln + nbreaks * (bytes_per_line - chars_per_line))
            res = re.sub(self.__removeWs, '', res)
        assert len(res) == ln, 'len(%s) != %d' % (res, ln)
        return res
