#!/usr/bin/env python
"""
spliced_read_recovery_performance.py

Outputs precision, recall, and related performance measurements of a spliced
aligner given a Flux-like BED file T with true spliced alignments and a SAM
file Y with the aligner's spliced alignments (read from stdin). Only reads that
overlap exon-exon junctions in the relevant (true) and retrieved
(aligner result) datasets are considered. Secondary alignments from the SAM
input are ignored.

If a coverage threshold is specified with -c, only reads that overlap at least
one exon-exon junction whose TRUE coverage is <= -c are considered.

THIS FILE DEPENDS ON DOOPLICITY AND RAIL; don't move it in the Rail repo.

Output (written to stdout)
----------------------------
Tab-delimited columns, where each row characterizes a distinct read alignment
    1. 1 if relevant else 0
    2. 1 if retrieved else 0
    3. comma-separated list of true coverages
    4. max true coverage
    5. min true coverage
    6. comma-separated list of retrieved coverages
    7. max retrieved coverage
    8. min retrieved coverage

(to stderr)
Two columns delimited by tabs, where the first column characterizes the numbers
in the second column. The first column is given below.

relevant instances
retrieved instances
intersection [of relevant and retrieved instances]
precision
recall
"""

import sys
import site
import os
import re
from collections import defaultdict

base_path = os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__))
                    )
utils_path = os.path.join(base_path, 'src', 'rna', 'utils')
src_path = os.path.join(base_path, 'src')

site.addsitedir(utils_path)
site.addsitedir(src_path)


def parsed_md(md):
    """ Divides an MD string up by boundaries between ^, letters, and numbers

        md: an MD string (example: 33A^CC).

        Return value: MD string split by boundaries described above.
    """
    md_to_parse = []
    md_group = [md[0]]
    for i, char in enumerate(md):
        if i == 0: continue
        if (re.match('[A-Za-z]', char) is not None) \
            != (re.match('[A-Za-z]', md[i-1]) is not None) or \
            (re.match('[0-9]', char) is not None) \
            != (re.match('[0-9]', md[i-1]) is not None):
            if md_group:
                md_to_parse.append(''.join(md_group))
            md_group = [char]
        else:
            md_group.append(char)
    if md_group:
        md_to_parse.append(''.join(md_group))
    return [char for char in md_to_parse if char != '0']


def indels_junctions_exons_mismatches(
        cigar, md, pos, seq, drop_deletions=False, junctions_only=False
):
    """ Finds indels, junctions, exons, mismatches from CIGAR, MD string, POS

        cigar: CIGAR string
        md: MD:Z string
        pos: position of first aligned base
        seq: read sequence
        drop_deletions: drops deletions from coverage vectors iff True
        junctions_only: does not populate mismatch list

        Return value: tuple
            (insertions, deletions, junctions, exons, mismatches).
        Insertions is a list of tuples (last genomic position before insertion,
                                 string of inserted bases). Deletions
            is a list of tuples (first genomic position of deletion,
                                 string of deleted bases). Junctions is a list
            of tuples (intron start position (inclusive),
                       intron end position (exclusive),
                       left_diplacement, right_displacement). Exons is a list
            of tuples (exon start position (inclusive),
                       exon end position (exclusive)). Mismatches is a list
            of tuples (genomic position of mismatch, read base)
    """
    insertions, deletions, junctions, exons, mismatches = [], [], [], [], []
    cigar = re.split(r'([MINDS])', cigar)[:-1]
    md = parsed_md(md)
    seq_size = len(seq)
    cigar_index, md_index, seq_index = 0, 0, 0
    max_cigar_index = len(cigar)
    while cigar_index != max_cigar_index:
        if cigar[cigar_index] == 0:
            cigar_index += 2
            continue
        if cigar[cigar_index + 1] == 'M':
            aligned_base_cap = int(cigar[cigar_index])
            aligned_bases = 0
            while True:
                try:
                    aligned_bases += int(md[md_index])
                    if aligned_bases <= aligned_base_cap:
                        md_index += 1
                except ValueError:
                    # Not an int, but should not have reached a deletion
                    assert md[md_index] != '^', '\n'.join(
                        ['cigar and md:',
                         ''.join(cigar), ''.join(md)]
                    )
                    if not junctions_only:
                        mismatches.append(
                            (pos + aligned_bases,
                             seq[seq_index + aligned_bases])
                        )
                    correction_length = len(md[md_index])
                    m_length = aligned_base_cap - aligned_bases
                    if correction_length > m_length:
                        md[md_index] = md[md_index][:m_length]
                        aligned_bases = aligned_base_cap
                    else:
                        aligned_bases += correction_length
                        md_index += 1
                if aligned_bases > aligned_base_cap:
                    md[md_index] = aligned_bases - aligned_base_cap
                    break
                elif aligned_bases == aligned_base_cap:
                    break
            # Add exon
            exons.append((pos, pos + aligned_base_cap))
            pos += aligned_base_cap
            seq_index += aligned_base_cap
        elif cigar[cigar_index + 1] == 'N':
            skip_increment = int(cigar[cigar_index])
            # Add junction
            junctions.append((pos, pos + skip_increment,
                              seq_index, seq_size - seq_index))
            # Skip region of reference
            pos += skip_increment
        elif cigar[cigar_index + 1] == 'I':
            # Insertion
            insert_size = int(cigar[cigar_index])
            insertions.append(
                (pos - 1, seq[seq_index:seq_index + insert_size])
            )
            seq_index += insert_size
        elif cigar[cigar_index + 1] == 'D':
            assert md[md_index] == '^', '\n'.join(
                ['cigar and md:',
                 ''.join(cigar), ''.join(md)]
            )
            # Deletion
            delete_size = int(cigar[cigar_index])
            md_delete_size = len(md[md_index + 1])
            assert md_delete_size >= delete_size
            deletions.append((pos, md[md_index + 1][:delete_size]))
            if not drop_deletions: exons.append((pos, pos + delete_size))
            if md_delete_size > delete_size:
                # Deletion contains a junction
                md[md_index + 1] = md[md_index + 1][delete_size:]
            else:
                md_index += 2
            # Skip deleted part of reference
            pos += delete_size
        else:
            # Soft clip
            assert cigar[cigar_index + 1] == 'S'
            # Advance seq_index
            seq_index += int(cigar[cigar_index])
        cigar_index += 2
    '''Merge exonic chunks/deletions; insertions/junctions could have chopped
    them up.'''
    new_exons = []
    last_exon = exons[0]
    for exon in exons[1:]:
        if exon[0] == last_exon[1]:
            # Merge ECs
            last_exon = (last_exon[0], exon[1])
        else:
            # Push last exon to new exon list
            new_exons.append(last_exon)
            last_exon = exon
    new_exons.append(last_exon)
    return insertions, deletions, junctions, new_exons, mismatches


class xstream(object):
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


def remove_temporary_directories(temp_dir_paths):
    """ Deletes temporary directory.

        temp_dir_paths: iterable of paths of temporary directories

        No return value.
    """
    for temp_dir_path in temp_dir_paths:
        try:
            shutil.rmtree(temp_dir_path,
                            ignore_errors=True)
        except Exception as e:
            # Don't know what's up, but forge on
            pass

def dummy_md_index(cigar):
    """ Creates dummy MD string from CIGAR in case of missing MD.

        cigar: cigar string

        Return value: dummy MD string
    """
    cigar = re.split(r'([MINDS])', cigar)[:-1]
    cigar_index = 0
    max_cigar_index = len(cigar)
    md = []
    while cigar_index != max_cigar_index:
        if cigar[cigar_index] == 0:
            cigar_index += 2
            continue
        if cigar[cigar_index+1] == 'M':
            try:
                if type(md[-1]) is int:
                    md[-1] += int(cigar[cigar_index])
                else:
                    md.append(int(cigar[cigar_index]))
            except IndexError:
                md.append(int(cigar[cigar_index]))
            cigar_index += 2
        elif cigar[cigar_index+1] in 'SIN':
            cigar_index += 2
        elif cigar[cigar_index+1] == 'D':
            md.extend(['^', 'A'*int(cigar[cigar_index])])
            cigar_index += 2
        else:
            raise RuntimeError(
                        'Accepted CIGAR characters are only in [MINDS].'
                    )
    return ''.join(str(el) for el in md)

"""
Here's what a few lines of the simulator BED file look like:

chr1	13362	13438	chr1:11869-14412W:ENST00000515242.2:2:1653:603:993/1	0	+	.	.	0,0,0	1	76	0
chr1	13677	13753	chr1:11869-14412W:ENST00000515242.2:2:1653:603:993/2	0	-	.	.	0,0,0	1	76	0
chr1	12693	13272	chr1:11869-14412W:ENST00000515242.2:3:1653:437:781/1	0	+	.	.	0,0,0	2	28,48	0,531
chr1	13465	13541	chr1:11869-14412W:ENST00000515242.2:3:1653:437:781/2	0	-	.	.	0,0,0	1	76	0
"""

def write_read_introns_from_bed_stream(bed_stream, output_stream,
                                        intron_counts,
                                        generous=False,
                                        instance=False):
    """ Writes output that maps QNAMES to exon-exon junctions overlapped.

        bed_stream: input stream containing lines of a BED file characterizing
            splice junctions.
        output_stream: where to write output. Each line takes the form:
            <read name><TAB><sorted list of intron starts and ends>
            <TAB>['t' for 'true']
        intron_counts: defaultdict(int) that counts number of simulated
            alignments overlapping exon-exon junction
        generous: True iff QNAMES should have the last two chars cut off

        No return value.
    """
    for line in bed_stream:
        introns = set()
        tokens = line.rstrip().split('\t')
        if len(tokens) < 12:
            continue  # BTL: what are these lines?
        chrom = tokens[0]
        chrom_start = int(tokens[1])
        chrom_end = int(tokens[2])
        name = tokens[3]  # this is a complicated string
        block_sizes = tokens[10].split(',')   # 1 and 2 above
        block_starts = tokens[11].split(',')  # 0 and 0,531 above
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
            intron_start = chrom_start + int(block_starts[i])     # BTL: so block starts are cumulative
            junctions.append(intron_start)                        # BTL: end of prev junction
            junctions.append(intron_start + int(block_sizes[i]))  # BTL: beginning of next junction
        # Final block characterizes junction on right side of intron
        junctions.append(chrom_start + int(block_starts[-1]))
        for i in xrange(len(junctions)/2):
            introns.add((junctions[2*i] + 1, junctions[2*i+1] + 1))
        if introns:
            introns = [chrom + ';'.join(['']
                                + [str(bound) for bound in intron])
                             for intron in sorted(list(introns))]
            if instance:
                print >>output_stream, '%s\t%s\tt' \
                    % (name[:-2] if generous else name, '\t'.join(introns))
            else:
                for intron in introns:
                    print >>output_stream, '%s;%s\t%s\tt' \
                    % (name[:-2] if generous else name,
                       intron,
                       intron)
        for intron in introns:
            intron_counts[intron] += 1

def write_read_introns_from_sam_stream(sam_stream, output_stream,
                                        retrieved_intron_counts,
                                        instance=False):
    """ Writes output that maps QNAMES to exon-exon junctions overlapped.
        BTL: also, populates the retrieved_intron_counts object.

        sam_stream: where to find retrieved alignments in SAM form
        output_stream: where to write output. Each line takes the form:
            <read name><TAB>RNAME<TAB><sorted list of intron starts and ends>
            <TAB>['r' for 'retrieved']
        retrieved_intron_counts: defaultdict(int) that counts number of
            retrieved alignments overlapping exon-exon junction

        No return value.
    """
    for line in sam_stream:
        if line[0] == '@': continue
        try:
            tokens = line.strip().split('\t')
            flag = int(tokens[1])
            if flag & 4:
                continue  # skip unmapped?  We won't miss true instances, will we?
            name = tokens[0]  # SAM: read name or chromosome name?
            rname = tokens[2]
            cigar = tokens[5]
            pos = int(tokens[3])
            seq = tokens[9]
            flag = int(tokens[1])
            if 'N' not in cigar or flag & 256:
                continue  # skip secondaries and ones with no splicing
            #md = [token[5:] for token in tokens if token[:5] == 'MD:Z:'][0]
            # BTL: it just wants introns back, so it doesn't matter that mismatches are getting ignored in the dummy md string
            _, _, introns, _, _ = indels_junctions_exons_mismatches(cigar,
                                        dummy_md_index(cigar), pos, seq)
            introns = [intron[:2] for intron in introns]  # just take first 2 fields of each
            introns = [rname 
                          + ';'.join([''] + [str(bound) for bound in intron])
                          for intron in sorted(list(introns))]
            if instance:
                for intron in introns:
                    retrieved_intron_counts[intron] += 1
                # Print read name, than a tab-separated list of introns covered,
                # each intron being a semicolon separated list of chromosome and
                # coordinates.  Then final column gets "r". (For retrieved?)
                print >>output_stream, '%s\t%s\tr' % (name, '\t'.join(introns))
            else:
                for intron in introns:
                    retrieved_intron_counts[intron] += 1
                    print >>output_stream, '%s;%s\t%s\tr' % (name,
                                                             intron,
                                                             intron)
        except IndexError:
            print >>sys.stderr, ('Error found on line: ' + line)
            raise

def print_instance(code, introns, intron_counts, retrieved_intron_counts):
    """ Prints information about intron coverage in an instance to stdout.

        An instance is a (read name, alignment) pair that appears in either or
        both of the original simulated data or the aligned data. A read that
        appears in the original data is relevant; a read that appears in the
        aligned data is retrieved. Here, "coverage" means the number of reads
        that overlap an exon-exon junction in the simulated or retrieved
        dataset.

        code: 'rel' if relevant but not retrieved
              'ret' if retrieved but not relevant
              'relret' if relevant and retrieved
        introns: list of introns crossed by instance
        intron_counts: maps elements of "introns" to their simulated coverages
        retrieved_intron_counts: maps elements of "introns" to their retrieved
            coverages

        Printed information is in the following format:
        1 if relevant else 0<TAB>1 if retrieved else 0<TAB>comma-separated
        list of coverages<TAB>max coverage<TAB>min coverage

        No return value.
    """
    if code == 'rel':
        relevant = True
        retrieved = False
    elif code == 'ret':
        relevant = False
        retrieved = True
    elif code == 'relret':
        relevant = True
        retrieved = True
    else:
        raise RuntimeError('Bad code: "%s"' % code)
    # BTL: does this properly hand 0 coverages?  Maybe it's a defaultdict?
    coverages = [intron_counts[intron] for intron in introns]
    retrieved_coverages = [retrieved_intron_counts[intron]
                            for intron in introns]
    print '%d\t%d\t%s\t%d\t%d\t%s\t%d\t%d' % (
            1 if relevant else 0,
            1 if retrieved else 0,
            ','.join([str(coverage) for coverage in coverages]),
            max(coverages),
            min(coverages),
            ','.join([str(retrieved_coverage) for retrieved_coverage
                        in retrieved_coverages]),
            max(retrieved_coverages),
            min(retrieved_coverages)
        )

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-t', '--true-introns-bed', type=str, required=True, 
        help='Full path of BED file containing true introns')
    parser.add_argument('-r', '--read-instance', action='store_const',
        const=True, default=False,
        help='Defines an instance as an alignment that overlaps exon-exon '
             'junctions rather than the event that an exon-exon junction '
             'is overlapped by an alignment')
    parser.add_argument('-g', '--generous', action='store_const', const=True,
        default=False,
        help='TopHat/STAR/HISAT cut off /1s and /2s from read names, even in '
             'unpaired mode. This loses information. In generous mode, '
             'this script provides extremely tight upper bounds on '
             'precision and recall for TopHat/STAR/HISAT')
    parser.add_argument('-c', '--coverage-threshold', type=int,
        required=False, default=None,
        help='Consider only reads that overlap at least '
             'one exon-exon junction whose TRUE coverage is <= -c')
    args = parser.parse_args(sys.argv[1:])
    import tempfile
    import atexit
    try:
        temp_dir_path = tempfile.mkdtemp(dir=os.path.dirname(
                                                    args.true_introns_bed
                                                )
                                            )
    except:
        temp_dir_path = tempfile.mkdtemp()
    #print >>sys.stderr, temp_dir_path
    atexit.register(remove_temporary_directories, [temp_dir_path])
    combined_file = os.path.join(temp_dir_path, 'combined.temp')
    # BTL: yes, these are defaultdicts
    intron_counts, retrieved_intron_counts = defaultdict(int), defaultdict(int)
    with open(combined_file, 'w') as combined_stream:
        with open(args.true_introns_bed) as true_introns_bed_stream:
            write_read_introns_from_bed_stream(true_introns_bed_stream,
                                                combined_stream,
                                                intron_counts,
                                                generous=args.generous,
                                                instance=args.read_instance)
        write_read_introns_from_sam_stream(sys.stdin, combined_stream,
                                            retrieved_intron_counts,
                                            instance=args.read_instance)
    import subprocess
    sorted_combined_file = os.path.join(temp_dir_path, 'combined.sorted.temp')
    subprocess.check_call(' '.join(['sort -T %s -k1,1' % temp_dir_path,
                                    combined_file, 
                                    '>', sorted_combined_file]),
                            bufsize=-1, shell=True)
    relevant = 0
    retrieved = 0
    relevant_and_retrieved = 0
    with open(sorted_combined_file) as sorted_combined_stream:
        for (name,), xpartition in xstream(sorted_combined_stream, 1):
            relevant_and_retrieved_instances = list(xpartition)
            ts = [instance[:-1] for instance 
                    in relevant_and_retrieved_instances
                    if instance[-1] == 't'
                    and (args.coverage_threshold is None
                         or any(
                    [intron_counts[intron] <= args.coverage_threshold
                     for intron in instance[:-1]]
                     ))]
            rs = [instance[:-1] for instance 
                    in relevant_and_retrieved_instances
                    if instance[-1] == 'r'
                    and (args.coverage_threshold is None
                         or any(
                    [intron_counts[intron] <= args.coverage_threshold
                     for intron in instance[:-1]]
                     ))]
            relevant_count = len(ts)
            retrieved_count = len(rs)
            assert relevant_count in [0, 1, 2]
            assert retrieved_count in [0, 1, 2]
            relevant += relevant_count
            retrieved += retrieved_count
            if retrieved_count == 2:
                '''For paired-end read output from TopHat and STAR, /1 and /2's
                are cut off from QNAMEs, so we must find the best assignment of
                retrieved alignments to relevant alignments.'''
                if relevant_count == 2:
                    aligned = [ts[0] == rs[0], ts[1] == rs[1]]
                    switched = [ts[0] == rs[1], rs[0] == ts[1]]
                    if aligned.count(True) >= switched.count(True):
                        matches = [(ts[0], rs[0]), (ts[1], rs[1])]
                    else:
                        matches = [(ts[0], rs[1]), (ts[1], rs[0])]
                    for match in matches:
                        if match[0] == match[1]:
                            relevant_and_retrieved += 1
                            print_instance('relret', match[0], intron_counts,
                                            retrieved_intron_counts)
                        else:
                            print_instance('rel', match[0], intron_counts,
                                            retrieved_intron_counts)
                            print_instance('ret', match[1], intron_counts,
                                            retrieved_intron_counts)
                elif relevant_count == 1:
                    try:
                        t_index = rs.index(ts[0])
                        for i, r in enumerate(rs):
                            if i == t_index:
                                print_instance('relret', r, intron_counts,
                                                retrieved_intron_counts)
                            else:
                                print_instance('ret', r, intron_counts,
                                                retrieved_intron_counts)
                        relevant_and_retrieved += 1
                    except ValueError:
                        # Relevant alignment not found
                        for r in rs:
                            print_instance('ret', r, intron_counts,
                                            retrieved_intron_counts)
                        print_instance('rel', ts[0], intron_counts,
                                        retrieved_intron_counts)
                else:
                    # relevant_count == 0
                    for r in rs:
                        print_instance('ret', r, intron_counts,
                                        retrieved_intron_counts)
            elif retrieved_count == 1:
                try:
                    r_index = ts.index(rs[0])
                    for i, t in enumerate(ts):
                        if i == r_index:
                            print_instance('relret', t, intron_counts,
                                            retrieved_intron_counts)
                        else:
                            print_instance('rel', t, intron_counts,
                                            retrieved_intron_counts)
                    relevant_and_retrieved += 1
                except ValueError:
                    # Retrieved alignment not found
                    for t in ts:
                        print_instance('rel', t, intron_counts,
                                        retrieved_intron_counts)
                    print_instance('ret', rs[0], intron_counts,
                                    retrieved_intron_counts)
            else:
                # No retrieved alignments
                for t in ts:
                    print_instance('rel', t, intron_counts,
                                    retrieved_intron_counts)
    precision = float(relevant_and_retrieved) / retrieved
    recall = float(relevant_and_retrieved) / relevant
    print >>sys.stderr, 'relevant instances\t%d' % relevant
    print >>sys.stderr, 'retrieved instances\t%d' % retrieved
    print >>sys.stderr, 'intersection\t%d' % relevant_and_retrieved
    print >>sys.stderr, 'precision\t%.12f' % precision
    print >>sys.stderr, 'recall\t%.12f' % recall
    print >>sys.stderr, ('fscore\t%.12f'
                            % (2 * precision * recall / (precision + recall)))
