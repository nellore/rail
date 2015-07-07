import os
import struct
import mmap
from operator import itemgetter
from collections import defaultdict
from bisect import bisect_right

class BowtieIndexReference(object):
    """
    Given prefix of a Bowtie index, parses the reference names, parses the
    extents of the unambiguous stretches, and memory-maps the file containing
    the unambiguous-stretch sequences.  get_stretch member function can
    retrieve stretches of characters from the reference, even if the stretch
    contains ambiguous characters.
    """

    def __init__(self, idx_prefix):

        # Open file handles
        if os.path.exists(idx_prefix + '.3.ebwt'):
            # Small index (32-bit offsets)
            fh1 = open(idx_prefix + '.1.ebwt', 'rb')  # for ref names
            fh3 = open(idx_prefix + '.3.ebwt', 'rb')  # for stretch extents
            fh4 = open(idx_prefix + '.4.ebwt', 'rb')  # for unambiguous sequence
            sz, struct_unsigned = 4, struct.Struct('I')
        else:
            raise RuntimeError('No Bowtie index files with prefix "%s"' % idx_prefix)

        #
        # Parse .1.bt2 file
        #
        one = struct.unpack('<i', fh1.read(4))[0]
        assert one == 1

        ln = struct_unsigned.unpack(fh1.read(sz))[0]
        line_rate = struct.unpack('<i', fh1.read(4))[0]
        lines_per_side = struct.unpack('<i', fh1.read(4))[0]
        _ = struct.unpack('<i', fh1.read(4))[0]
        ftab_chars = struct.unpack('<i', fh1.read(4))[0]
        _ = struct.unpack('<i', fh1.read(4))[0]

        nref = struct_unsigned.unpack(fh1.read(sz))[0]
        # get ref lengths
        reference_length_list = []
        for i in xrange(nref):
            reference_length_list.append(struct.unpack('<i', fh1.read(sz))[0])

        nfrag = struct_unsigned.unpack(fh1.read(sz))[0]
        # skip rstarts
        fh1.seek(nfrag * sz * 3, 1)

        # skip ebwt
        bwt_sz = ln // 4 + 1
        line_sz = 1 << line_rate
        side_sz = line_sz * lines_per_side
        side_bwt_sz = side_sz - 8
        num_side_pairs = (bwt_sz + (2*side_bwt_sz) - 1) // (2*side_bwt_sz)
        ebwt_tot_len = num_side_pairs * 2 * side_sz
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

        running_unambig, running_length = 0, 0
        self.recs = defaultdict(list)
        self.offset_in_ref = defaultdict(list)
        self.unambig_preceding = defaultdict(list)
        length = {}

        ref_id, ref_namenrecs_added = 0, None
        for i in xrange(nrecs):
            off = struct_unsigned.unpack(fh3.read(sz))[0]
            ln = struct_unsigned.unpack(fh3.read(sz))[0]
            first_of_chromosome = ord(fh3.read(1)) != 0
            if first_of_chromosome:
                if i > 0:
                    length[ref_name] = running_length
                ref_name = refnames[ref_id]
                ref_id += 1
                running_length = 0
            assert ref_name is not None
            self.recs[ref_name].append((off, ln, first_of_chromosome))
            self.offset_in_ref[ref_name].append(running_length)
            self.unambig_preceding[ref_name].append(running_unambig)
            running_length += (off + ln)
            running_unambig += ln

        length[ref_name] = running_length
        assert nrecs == sum(map(len, self.recs.itervalues()))

        #
        # Memory-map the .4.bt2 file
        #
        ln_bytes = (running_unambig + 3) // 4
        self.fh4mm = mmap.mmap(fh4.fileno(), ln_bytes, flags=mmap.MAP_SHARED, prot=mmap.PROT_READ)

        # These are per-reference
        self.length = length
        self.refnames = refnames

        # To facilitate sorting reference names in order of descending length
        sorted_rnames = sorted(self.length.items(),
                               key=lambda x: itemgetter(1)(x), reverse=True)
        '''A case-sensitive sort is also necessary here because new versions of
        bedGraphToBigWig complain on encountering a nonlexicographic sort
        order.'''
        lexicographically_sorted_rnames = sorted(self.length.items(),
                                                    key=lambda x:
                                                    itemgetter(0)(x))
        self.rname_to_string, self.l_rname_to_string = {}, {}
        self.string_to_rname, self.l_string_to_rname = {}, {}
        for i, (rname, _) in enumerate(sorted_rnames):
            rname_string = ('%012d' % i)
            self.rname_to_string[rname] = rname_string
            self.string_to_rname[rname_string] = rname
        for i, (rname, _) in enumerate(lexicographically_sorted_rnames):
            rname_string = ('%012d' % i)
            self.l_rname_to_string[rname] = rname_string
            self.l_string_to_rname[rname_string] = rname
        # Handle unmapped reads
        unmapped_string = ('%012d' % len(sorted_rnames))
        self.rname_to_string['*'] = unmapped_string
        self.string_to_rname[unmapped_string] = '*'

        # For compatibility
        self.rname_lengths = self.length

    def get_stretch(self, ref_id, ref_off, count):
        """
        Return a stretch of characters from the reference, retrieved
        from the Bowtie index.

        @param ref_id: name of ref seq, up to & excluding whitespace
        @param ref_off: offset into reference, 0-based
        @param count: # of characters
        @return: string extracted from reference
        """
        assert ref_id in self.recs
        # Account for negative reference offsets by padding with Ns
        N_count = min(abs(min(ref_off, 0)), count)
        stretch = ['N'] * N_count
        count -= N_count
        if not count: return ''.join(stretch)
        ref_off = max(ref_off, 0)
        starting_rec = bisect_right(self.offset_in_ref[ref_id], ref_off) - 1
        assert starting_rec >= 0
        off = self.offset_in_ref[ref_id][starting_rec]
        buf_off = self.unambig_preceding[ref_id][starting_rec]
        # Naive to scan these records linearly; obvious speedup is binary search
        for rec in self.recs[ref_id][starting_rec:]:
            off += rec[0]
            while ref_off < off and count > 0:
                stretch.append('N')
                count -= 1
                ref_off += 1
            if count == 0:
                break
            if ref_off < off + rec[1]:
                # stretch extends through part of the unambiguous stretch
                buf_off += (ref_off - off)
            else:
                buf_off += rec[1]
            off += rec[1]
            while ref_off < off and count > 0:
                buf_elt = buf_off >> 2
                shift_amt = (buf_off & 3) << 1
                stretch.append(
                    'ACGT'[(ord(self.fh4mm[buf_elt]) >> shift_amt) & 3]
                )
                buf_off += 1
                count -= 1
                ref_off += 1
            if count == 0:
                break
        # If the requested stretch went past the last unambiguous
        # character in the chromosome, pad with Ns
        while count > 0:
            count -= 1
            stretch.append('N')
        return ''.join(stretch)


def which(program):
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

        class TestBowtieIndexReference(unittest.TestCase):

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
                assert which('bowtie-build') is not None
                os.system('bowtie-build %s %s >/dev/null' % (self.fa_fn_1, self.fa_fn_1))

            def tearDown(self):
                rmtree(self.tmpdir)

            def test1(self):
                ref = BowtieIndexReference(self.fa_fn_1)
                self.assertEqual('ACGTACGTAC', ref.get_stretch('short_name1', 0, 10))
                self.assertEqual('ACGTACGTAC', ref.get_stretch('short_name1', 40, 10))
                self.assertEqual('ANNNNNNNNN', ref.get_stretch('short_name1', 80, 10))

                self.assertEqual('CAGTCAGTCA', ref.get_stretch('short_name2', 0, 10))
                self.assertEqual('CAGTCAGTCA', ref.get_stretch('short_name2', 40, 10))
                self.assertEqual('NNNNNNNNNN', ref.get_stretch('short_name2', 80, 10))

                self.assertEqual('CANNNNNNNN', ref.get_stretch('short_name3', 0, 10))

            def test2(self):
                ref = BowtieIndexReference(self.fa_fn_1)
                self.assertEqual('CAGTCAGTCA', ref.get_stretch('short_name2', 0, 10))
                self.assertEqual('AGTCAGTCAGT', ref.get_stretch('short_name2', 1, 11))
                self.assertEqual('GTCAGTCAGTCA', ref.get_stretch('short_name2', 2, 12))
                self.assertEqual('TCAGTCAGTCAGT', ref.get_stretch('short_name2', 3, 13))

                self.assertEqual('TACGTACGTA', ref.get_stretch('short_name1', 71, 10))
                self.assertEqual('ACGTACGTANN', ref.get_stretch('short_name1', 72, 11))
                self.assertEqual('CGTACGTANNNN', ref.get_stretch('short_name1', 73, 12))

            def test3(self):
                ref = BowtieIndexReference(self.fa_fn_1)
                self.assertEqual('NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN', ref.get_stretch('short_name4', 0, 40))
                self.assertEqual('NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNA', ref.get_stretch('short_name4', 1, 40))
                self.assertEqual('AAAA', ref.get_stretch('short_name4', 41, 4))
                self.assertEqual('NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTT', ref.get_stretch('short_name4', 240, 42))

            def test4(self):
                ref = BowtieIndexReference(self.fa_fn_1)
                # Test that all refname lengths are accurate
                self.assertEqual(ref.length['short_name1'], 81)
                self.assertEqual(ref.length['short_name4'], 282)
                self.assertEqual(ref.length['short_name2'], 80)
                self.assertEqual(ref.length['short_name3'], 2)

            def test_off_reference_values(self):
                ref = BowtieIndexReference(self.fa_fn_1)
                self.assertEqual('NNNACG', ref.get_stretch('short_name1', -3, 6))
                self.assertEqual('NNNNN', ref.get_stretch('short_name1', -20, 5))
                self.assertEqual('NNNNNNNNN', ref.get_stretch('short_name1', 85, 9))
                self.assertEqual('ANNNNNNNN', ref.get_stretch('short_name1', 80, 9))

        unittest.main(argv=[sys.argv[0]])
