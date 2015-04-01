"""
count_introns.py

Uses bowtie2-inspect to get RNAMEs from transcript fragment index written by
 Rail-RNA and count the number of introns in it. Specification of RNAMEs from
this index is in intron_fasta.py.

Once introns are grabbed, counts how many have canonical and noncanonical
splice sites using BowtieIndexReference, a class from Rail-RNA pasted here so
the evaluation code is independent of Rail-RNA's source.

This script may be useful to users who want to interpret the results in
the transcript_index subdirectory of the output directory of a Rail-RNA run;
it provides a list of the introns uncovered by Rail-RNA before any realignment
steps.
"""

import os
import struct
import mmap
from operator import itemgetter
from bisect import bisect_right
from collections import defaultdict

def introns_from_bed(bed):
    """ Converts BED to dictionary that maps RNAMES to sets of introns.

        bed: input BED file characterizing splice junctions

        Return value: a dictionary. Each key is an RNAME, typically a
            chromosome, and its corresponding value is a set of tuples, each
            denoting an intron on RNAME. Each tuple is of the form
            (start position, end position).
    """
    introns = set()
    with open(bed) as bed_stream:
        for line in bed_stream:
            tokens = line.rstrip().split('\t')
            if len(tokens) < 12:
                continue
            chrom = tokens[0]
            chrom_start = int(tokens[1])
            chrom_end = int(tokens[2])
            block_sizes = tokens[10].split(',')
            block_starts = tokens[11].split(',')
            # Handle trailing commas
            try:
                int(block_sizes[-1])
            except ValueError:
                block_sizes = block_sizes[:-1]
            try:
                int(block_starts[-1])
            except ValueError:
                block_starts = block_starts[:-1]
            block_count = len(block_sizes)
            if block_count < 2:
                # No introns
                continue
            assert block_count == len(block_starts)
            junctions = []
            # First block characterizes junction on left side of intron
            junctions.append(chrom_start + int(block_starts[0]) 
                                    + int(block_sizes[0]))
            for i in xrange(1, block_count - 1):
                # Any intervening blocks characterize two junctions
                intron_start = chrom_start + int(block_starts[i])
                junctions.append(intron_start)
                junctions.append(intron_start + int(block_sizes[i]))
            # Final block characterizes junction on right side of intron
            junctions.append(chrom_start + int(block_starts[-1]))
            for i in xrange(len(junctions)/2):
                introns.add((chrom, junctions[2*i]+1, junctions[2*i+1]+1))
    return introns

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
            fh4 = open(idx_prefix + '.4.ebwt', 'rb')  # for unambiguous seq
            sz, struct_unsigned = 4, struct.Struct('I')
        else:
            raise RuntimeError('No Bowtie index files with prefix "%s"'
                                    % idx_prefix)

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
        self.fh4mm = mmap.mmap(fh4.fileno(), ln_bytes,
                                flags=mmap.MAP_SHARED,
                                prot=mmap.PROT_READ)

        # These are per-reference
        self.length = length
        self.refnames = refnames

        # To facilitate sorting reference names in order of descending length
        sorted_rnames = sorted(self.length.items(),
                               key=lambda x: itemgetter(1)(x), reverse=True)
        self.rname_to_string = {}
        self.string_to_rname = {}
        for i, (rname, _) in enumerate(sorted_rnames):
            rname_string = ('%012d' % i)
            self.rname_to_string[rname] = rname_string
            self.string_to_rname[rname_string] = rname
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
        stretch = []
        starting_rec = bisect_right(self.offset_in_ref[ref_id], ref_off) - 1
        assert starting_rec >= 0
        off = self.offset_in_ref[ref_id][starting_rec]
        buf_off = self.unambig_preceding[ref_id][starting_rec]
        '''Naive to scan these records linearly; obvious speedup is binary
        search'''
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

if __name__ == '__main__':
    import argparse
    import subprocess

    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--bowtie2-inspect', type=str, required=False,
            default='bowtie2-inspect',
            help='Path to Bowtie 2 executable'
        )
    parser.add_argument('--bowtie1-idx', type=str, required=True,
            help=('Path to basename of Bowtie 1 index with genome to which '
                  'FASTQs were aligned')
        )
    parser.add_argument('--bowtie2-idx', type=str, required=True,
            help=('Path to basename of Bowtie 2 index containing transcript '
                  'fragments')
        )
    parser.add_argument('-t', '--true-introns-bed-dir', type=str,
        required=False,
        default=None,
        help=('Path to directory containing Flux Simulator BEDs. The only '
              'BEDs in the directory (subdirectories excluded) must be '
              'from the Flux sim of 112 bioreps. Exclude if precision/recall '
              'are not to be computed.'))
    parser.add_argument('--print-introns-to-stderr', action='store_const',
                            const=True, default=False,
                            help=('Prints introns to stderr in format '
                                  'RNAME <TAB> sense <TAB> intron start <TAB> '
                                  'intron end. Coordinates are '
                                  '1-based; start coordinate is inclusive '
                                  'while end coordinate is exclusive.'))

    args = parser.parse_args()

    inspect_process = subprocess.Popen(
                            [args.bowtie2_inspect, '-n', args.bowtie2_idx],
                            stdout=subprocess.PIPE
                        )
    introns = set()
    for line in inspect_process.stdout:
        rname_and_sense, seq_start, subseq_sizes, intron_sizes, _, _ = \
            line.split('\x1d')
        seq_start = int(seq_start)
        subseq_sizes = [int(size) for size in subseq_sizes.split(',')]
        intron_sizes = [int(size) for size in intron_sizes.split(',')]
        assert len(subseq_sizes) == len(intron_sizes) + 1
        start = seq_start + subseq_sizes[0]
        for i, size in enumerate(subseq_sizes[1:]):
            introns.add((
                       rname_and_sense[:-1],
                       start, start + intron_sizes[i]
                    )
            )
            start += (size + intron_sizes[i])
    inspect_process.wait()
    possible_combos = set([('GT', 'AG'), ('CT', 'AC'),
                            ('GC', 'AG'), ('CT', 'GC'),
                            ('AT', 'AC'), ('GT', 'AT')])
    canonical = set([('GT', 'AG'), ('CT', 'AC')])
    less_canonical = set([('GC', 'AG'), ('CT', 'GC')])
    much_less_canonical = set([('AT', 'AC'), ('GT', 'AT')])
    reference_index = BowtieIndexReference(args.bowtie1_idx)

    canonicals, less_canonicals, much_less_canonicals = 0, 0, 0
    for rname, start, end in introns:
        left = reference_index.get_stretch(rname, start - 1, 2)
        right = reference_index.get_stretch(rname, end - 3, 2)
        assert (left, right) in possible_combos, (left, right)
        if (left, right) in canonical:
            canonicals += 1
        elif (left, right) in less_canonical:
            less_canonicals += 1
        elif (left, right) in much_less_canonical:
            much_less_canonicals += 1

    intron_count = len(introns)
    print 'GT-AG count: %d\tproportion: %08f' % (canonicals, float(canonicals)
                                                                / intron_count)
    print 'GC-AG count: %d\tproportion: %08f' % (less_canonicals,
                                        float(less_canonicals) / intron_count)
    print 'AT-AC count: %d\tproportion: %08f' % (much_less_canonicals,
                                        float(much_less_canonicals)
                                                                / intron_count)
    if args.true_introns_bed_dir is not None:
        # Read Flux BEDs
        true_introns = set()
        def add_sets(list_of_sets):
            """ For updating set with sets in list

                list_of_sets: list to sets to add to true_introns

                No return value.
            """
            global true_introns
            for item in list_of_sets:
                true_introns.update(item)
        import glob
        import multiprocessing
        pool = multiprocessing.Pool(multiprocessing.cpu_count() - 1)
        pool.map_async(
                    introns_from_bed,
                    glob.glob(
                            os.path.join(args.true_introns_bed_dir, '*.bed')
                        ),
                    callback=add_sets
                )
        pool.close()
        pool.join()
        retrieved = intron_count
        relevant = len(true_introns)
        relevant_and_retrieved = len(introns.intersection(true_introns))
        print 'true intron count\t%d' % relevant
        print 'retrieved intron count\t%d' % retrieved
        print 'overlap\t%d' % relevant_and_retrieved
        print 'precision\t%.9f' % (float(relevant_and_retrieved) / retrieved)
        print 'recall\t%.9f' % (float(relevant_and_retrieved) / relevant)
    else:
        print 'intron count\t%d' % intron_count
    if args.print_introns_to_stderr:
        introns = sorted(
                            list(introns),
                            key=lambda intron: (intron[0][:-1], intron[2])
                    )
        for rname, start, end in introns:
            print '\t'.join([rname[:-1], str(start), str(end)])