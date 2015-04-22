#!/usr/bin/env python
"""
Rail-RNA-bed_pre
Follows Rail-RNA-realign_reads
Precedes Rail-RNA-bed / Rail-RNA-tsv

Reducer for MapReduce pipelines that sums coverage of indels/introns and
computes other statistics by sample from Rail-RNA-realign.

Input (read from stdin)
----------------------------
Tab-delimited input tuple columns:
1. 'I', 'D', or 'N' for insertion, deletion, or intron line
2. Number string representing RNAME
3. Start position (Last base before insertion, first base of deletion,
                    or first base of intron)
4. End position (Last base before insertion, last base of deletion (exclusive),
                    or last base of intron (exclusive))
5. '+' or '-' indicating which strand is the sense strand for introns,
   inserted sequence for insertions, or deleted sequence for deletions
6. Sample index
----Next fields are for introns only; they are '\x1c' for indels----
7. Number of nucleotides between 5' end of intron and 5' end of read from which
    it was inferred, ASSUMING THE SENSE STRAND IS THE FORWARD STRAND. That is,
    if the sense strand is the reverse strand, this is the distance between the
    3' end of the read and the 3' end of the intron.
8. Number of nucleotides between 3' end of intron and 3' end of read from which
    it was inferred, ASSUMING THE SENSE STRAND IS THE FORWARD STRAND.
--------------------------------------------------------------------
9. Number of instances of intron, insertion, or deletion in sample; this is
    always +1 before bed_pre combiner/reducer

Input is partitioned by fields 1-5 and sorted by field 6.

Hadoop output (written to stdout)
----------------------------
Tab-delimited output tuple columns (bed):
1. 'I', 'D', or 'N' for insertion, deletion, or intron line
2. Sample index
3. Number string representing RNAME (+ '+ or -' if intron; same as field 6)
4. Start position (Last base before insertion, first base of deletion,
                    or first base of intron)
5. End position (Last base before insertion, last base of deletion (exclusive),
                    or last base of intron (exclusive))
6. '+' or '-' indicating which strand is the sense strand for introns,
   inserted sequence for insertions, or deleted sequence for deletions
----Next fields are for introns only; they are '\x1c' for indels----
7. MAX number of nucleotides between 5' end of intron and 5' end of read from
    which it was inferred, ASSUMING THE SENSE STRAND IS THE FORWARD STRAND.
    That is, if the sense strand is the reverse strand, this is the distance
    between the 3' end of the read and the 3' end of the intron.
8. MAX number of nucleotides between 3' end of intron and 3' end of read from
    which it was inferred, ASSUMING THE SENSE STRAND IS THE FORWARD STRAND.
9. MAXIMIN (number of nucleotides between 5' end of intron and 5' end of read,
            number of nucleotides between 3' end of intron and 3' end of read);
   min is between the args above; max is across reads.

Tab-delimited output tuple columns (collect)
1. '0' if insertion, '1' if deletion, or '2' if intron line
2. Number string representing RNAME (+ '+ or -' if intron; same as field 6)
3. Start position (Last base before insertion, first base of deletion,
                    or first base of intron)
4. End position (Last base before insertion, last base of deletion (exclusive),
                    or last base of intron (exclusive))
5. '+' or '-' indicating which strand is the sense strand for introns,
   inserted sequence for insertions, or deleted sequence for deletions
6. Coverage of feature for sample with index N
...
N + 6. Coverage of feature in sample with index N
--------------------------------------------------------------------
10. SUMMED number of instances of intron, insertion, or deletion in sample

OUTPUT COORDINATES ARE 1-INDEXED.
"""
import os
import sys
import argparse
import site

base_path = os.path.abspath(
                    os.path.dirname(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__)))
                    )
                )
utils_path = os.path.join(base_path, 'rna', 'utils')
site.addsitedir(utils_path)
site.addsitedir(base_path)
from dooplicity.tools import xstream

# Print file's docstring if -h is invoked
parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
args = parser.parse_args()

import time
import itertools
start_time = time.time()
input_line_count, output_line_count = 0, 0

for (line_type, rname, pos, end_pos, strand_or_seq), xpartition in xstream(
            sys.stdin, 5
        ):
    collect_line = [rname, pos, end_pos, strand_or_seq]
    i = 0
    if line_type == 'N':
        for sample_index, data in itertools.groupby(xpartition, 
                                                    key=lambda val: val[0]):
            sample_index = int(sample_index)
            while i != sample_index:
                # Write 0 coverage for sample indexes reporting 0 introns
                collect_line.append('0')
                i += 1
            coverage_sum = 0
            max_left_displacement, max_right_displacement = None, None
            maximin_displacement = None
            for _, left_displacement, right_displacement, coverage in data:
                input_line_count += 1
                left_displacement = int(left_displacement)
                right_displacement = int(right_displacement)
                max_left_displacement = max(left_displacement,
                                            max_left_displacement)
                max_right_displacement = max(right_displacement,
                                             max_right_displacement)
                maximin_displacement = max(
                        min(left_displacement, right_displacement),
                        maximin_displacement
                    )
                coverage_sum += int(coverage)
            assert max_left_displacement is not None
            assert max_right_displacement is not None
            assert maximin_displacement is not None
            print 'bed\tN\t%d\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d' % (
                    sample_index, rname, pos, end_pos, strand_or_seq,
                    max_left_displacement, max_right_displacement,
                    maximin_displacement, coverage_sum
                )
            collect_line.append(str(coverage_sum))
            i += 1
            output_line_count += 1
        sys.stdout.write('collect\t2\t')
    else:
        assert line_type in 'ID'
        for sample_index, data in itertools.groupby(xpartition, 
                                                    key=lambda val: val[0]):
            sample_index = int(sample_index)
            while i != sample_index:
                # Write 0 coverage for sample indexes reporting 0 introns
                collect_line.append('0')
                i += 1
            coverage_sum = 0
            for _, _, _, coverage in data:
                input_line_count += 1
                coverage_sum += int(coverage)
            print  'bed\t%s\t%s\t%s\t%s\t%s\t%s\t\x1c\t\x1c\t\x1c\t%d' % (
                    line_type, sample_index, rname, pos, end_pos,
                    strand_or_seq, coverage_sum
                )
            collect_line.append(str(coverage_sum))
            i += 1
            output_line_count += 1
        if line_type == 'I':
            sys.stdout.write('collect\t0\t')
        else:
            sys.stdout.write('collect\t1\t')
    print '\t'.join(collect_line)

print >>sys.stderr, 'DONE with bed_pre.py; in/out =%d/%d; time=%0.3f s' \
                        % (input_line_count, output_line_count,
                            time.time() - start_time)
