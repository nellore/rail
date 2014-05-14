import os
import struct
import mmap


class Bowtie2IndexReference(object):
    """
    Given prefix of a Bowtie 2 index, parses the reference names, parses the
    extents of the unambiguous stretches, and memory-maps the file containing
    the unambiguous-stretch sequences.  get_stretch member function can
    retrieve stretches of characters from the reference, even if the stretch
    contains ambiguous characters.
    """

    def __init__(self, idx_prefix):

        # Open file handles
        if os.path.exists(idx_prefix + '.3.bt2'):
            # Small index (32-bit offsets)
            fh1 = open(idx_prefix + '.1.bt2', 'rb')  # for ref names
            fh3 = open(idx_prefix + '.3.bt2', 'rb')  # for stretch extents
            fh4 = open(idx_prefix + '.4.bt2', 'rb')  # for unambiguous sequence
            sz, struct_unsigned = 4, struct.Struct('I')
        elif os.path.exists(idx_prefix + '.3.bt2l'):
            # Large index (64-bit offsets)
            fh1 = open(idx_prefix + '.1.bt2l', 'rb')  # for ref names
            fh3 = open(idx_prefix + '.3.bt2l', 'rb')  # for stretch extents
            fh4 = open(idx_prefix + '.4.bt2l', 'rb')  # for unambiguous sequence
            sz, struct_unsigned = 8, struct.Struct('Q')
        else:
            raise RuntimeError('No Bowtie 2 index files with prefix "%s"' % idx_prefix)

        #
        # Parse .1.bt2 file
        #
        one = struct.unpack('<i', fh1.read(4))[0]
        assert one == 1

        ln = struct_unsigned.unpack(fh1.read(sz))[0]
        line_rate = struct.unpack('<i', fh1.read(4))[0]
        _ = struct.unpack('<i', fh1.read(4))[0]
        _ = struct.unpack('<i', fh1.read(4))[0]
        ftab_chars = struct.unpack('<i', fh1.read(4))[0]
        _ = struct.unpack('<i', fh1.read(4))[0]

        nref = struct_unsigned.unpack(fh1.read(sz))[0]
        # skip ref lengths
        fh1.seek(nref * sz, 1)

        nfrag = struct_unsigned.unpack(fh1.read(sz))[0]
        # skip rstarts
        fh1.seek(nfrag * sz * 3, 1)

        # skip ebwt
        bwt_sz = ln // 4 + 1
        line_sz = 1 << line_rate
        side_sz = line_sz
        side_bwt_sz = side_sz - (sz * 4)
        num_sides = (bwt_sz + side_bwt_sz - 1) // side_bwt_sz
        ebwt_tot_len = num_sides * side_sz
        fh1.seek(ebwt_tot_len, 1)

        # skip zOff
        fh1.seek(sz, 1)

        # skip fchr
        fh1.seek(5 * sz, 1)

        # skip ftab
        ftab_len = (1 << (ftab_chars * 2)) + 1
        fh1.seek(ftab_len * sz, 1)

        # skip eftab
        eftab_len = ftab_chars * 2
        fh1.seek(eftab_len * sz, 1)

        refnames = []
        while True:
            refname = fh1.readline()
            if len(refname) == 0 or ord(refname[0]) == 0:
                break
            refnames.append(refname.split()[0])
        assert len(refnames) == nref

        #
        # Parse .3.bt2 file
        #
        one = struct.unpack('<i', fh3.read(4))[0]
        assert one == 1

        nrecs = struct_unsigned.unpack(fh3.read(sz))[0]

        running_unambig_preceding, running_length = 0, 0
        recs = []
        starting_offsets = []
        unambig_preceding, lengths = [], []

        for i in xrange(nrecs):
            off = struct_unsigned.unpack(fh3.read(sz))[0]
            ln = struct_unsigned.unpack(fh3.read(sz))[0]
            first = ord(fh3.read(1)) != 0
            recs.append((off, ln, first))
            if first:
                starting_offsets.append(len(recs)-1)
                unambig_preceding.append(running_unambig_preceding)
                if i > 0:
                    lengths.append(running_length)
                running_length = 0
            running_length += (off + ln)
            running_unambig_preceding += ln

        lengths.append(running_length)
        starting_offsets.append(len(recs))
        tot_unambig_len = running_unambig_preceding
        assert len(recs) == nrecs

        #
        # Memory-map the .4.bt2 file
        #
        ln_bytes = (tot_unambig_len + 3) // 4
        self.fh4mm = mmap.mmap(fh4.fileno(), ln_bytes, flags=mmap.MAP_SHARED, prot=mmap.PROT_READ)

        # These are per-unambiguous-stretch
        self.recs = recs

        # These are per-reference
        self.unambig_preceding = unambig_preceding
        self.lengths = lengths
        self.starting_offsets = starting_offsets
        self.refnames = refnames
        self.ref_id_to_offset = {self.refnames[i]: i for i in xrange(len(self.refnames))}

    def get_stretch(self, ref_id, ref_off, count):
        assert ref_id in self.ref_id_to_offset
        ref_idx = self.ref_id_to_offset[ref_id]
        rec_i, rec_f = self.starting_offsets[ref_idx], self.starting_offsets[ref_idx + 1]
        cur, off = 0, 0
        buf_off = self.unambig_preceding[ref_idx]
        stretch = []
        # Naive to scan these records linearly; obvious speedup is binary search
        for i in xrange(rec_i, rec_f):
            off += self.recs[i][0]
            while ref_off < off and count > 0:
                stretch.append('N')
                count -= 1
                ref_off += 1
            if count == 0:
                break
            if ref_off < off + self.recs[i][1]:
                # stretch extends through part of the unambiguous stretch
                buf_off += (ref_off - off)
            else:
                buf_off += self.recs[i][1]
            off += self.recs[i][1]
            while ref_off < off and count > 0:
                buf_elt = buf_off >> 2
                shift_amt = (buf_off & 3) << 1
                stretch.append('ACGT'[(ord(self.fh4mm[buf_elt]) >> shift_amt) & 3])
                buf_off += 1
                count -= 1
                ref_off += 1
            if count == 0:
                break
        while count > 0:
            count -= 1
            stretch.append('N')
        return ''.join(stretch)


def which(program):
    import os

    def is_exe(fp):
        return os.path.isfile(fp) and os.access(fp, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

if __name__ == "__main__":

    import sys
    import unittest
    import argparse
    from tempfile import mkdtemp
    from shutil import rmtree

    parser = argparse.ArgumentParser()
    parser.add_argument('--test', action='store_const', const=True, default=False, help='Do unit tests')

    args = parser.parse_args()

    if args.test:
        import unittest

        class TestBowtie2IndexReference(unittest.TestCase):

            def setUp(self):
                self.tmpdir = mkdtemp()
                self.fa_fn_1 = os.path.join(self.tmpdir, 'tmp1.fa')
                with open(self.fa_fn_1, 'w') as fh:
                    fh.write('''>short_name1 with some stuff after whitespace
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
A
>short_name4
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
TT
>short_name2 with some stuff after whitespace
CAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGT
CAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGT
>short_name3 with some stuff after whitespace
CA
''')
                assert which('bowtie2-build') is not None
                os.system('bowtie2-build %s %s >/dev/null 2>/dev/null' % (self.fa_fn_1, self.fa_fn_1))

            def tearDown(self):
                rmtree(self.tmpdir)

            def test1(self):
                ref = Bowtie2IndexReference(self.fa_fn_1)
                self.assertEqual('ACGTACGTAC', ref.get_stretch('short_name1', 0, 10))
                self.assertEqual('ACGTACGTAC', ref.get_stretch('short_name1', 40, 10))
                self.assertEqual('ANNNNNNNNN', ref.get_stretch('short_name1', 80, 10))

                self.assertEqual('CAGTCAGTCA', ref.get_stretch('short_name2', 0, 10))
                self.assertEqual('CAGTCAGTCA', ref.get_stretch('short_name2', 40, 10))
                self.assertEqual('NNNNNNNNNN', ref.get_stretch('short_name2', 80, 10))

                self.assertEqual('CANNNNNNNN', ref.get_stretch('short_name3', 0, 10))

            def test2(self):
                ref = Bowtie2IndexReference(self.fa_fn_1)
                self.assertEqual('CAGTCAGTCA', ref.get_stretch('short_name2', 0, 10))
                self.assertEqual('AGTCAGTCAGT', ref.get_stretch('short_name2', 1, 11))
                self.assertEqual('GTCAGTCAGTCA', ref.get_stretch('short_name2', 2, 12))
                self.assertEqual('TCAGTCAGTCAGT', ref.get_stretch('short_name2', 3, 13))

                self.assertEqual('TACGTACGTA', ref.get_stretch('short_name1', 71, 10))
                self.assertEqual('ACGTACGTANN', ref.get_stretch('short_name1', 72, 11))
                self.assertEqual('CGTACGTANNNN', ref.get_stretch('short_name1', 73, 12))

            def test3(self):
                ref = Bowtie2IndexReference(self.fa_fn_1)
                self.assertEqual('NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN', ref.get_stretch('short_name4', 0, 40))
                self.assertEqual('NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNA', ref.get_stretch('short_name4', 1, 40))
                self.assertEqual('AAAA', ref.get_stretch('short_name4', 41, 4))
                self.assertEqual('NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTT', ref.get_stretch('short_name4', 240, 42))

        unittest.main(argv=[sys.argv[0]])
