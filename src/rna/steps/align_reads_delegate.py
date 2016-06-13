"""
align_reads_delegate.py 

Output of Bowtie 2 from align_reads.py is streamed to this script to obtain
first and second-pass output. See align_reads.py for output format information.
"""

import sys
import os
import site
import string
from collections import defaultdict
import re

if '--test' in sys.argv:
    print("No unit tests")
    #unittest.main(argv=[sys.argv[0]])
    sys.exit(0)

base_path = os.path.abspath(
                    os.path.dirname(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__)))
                    )
                )
utils_path = os.path.join(base_path, 'rna', 'utils')
site.addsitedir(utils_path)
site.addsitedir(base_path)

from dooplicity.tools import xstream, xopen
from alignment_handlers import AlignmentPrinter
import bowtie_index
import bowtie
import manifest
import partition
import group_reads
from encode import decode_sequence

_reversed_complement_translation_table = string.maketrans('ATCG', 'TAGC')
_output_line_count = 0
_polyA = set(['A'])

def print_readletized_output(seq, sample_indexes,
        reversed_complement_sample_indexes, seq_id, cap_sizes,
        output_stream=sys.stdout, min_readlet_size=8, max_readlet_size=25,
        readlet_interval=5, verbose=False, no_polyA=False):
    """ Readletizes the unique read sequence seq.

        The function divides a sequence seq into several overlapping segments,
        or readlets. These readlets will later be aligned and used to infer the
        junctions that original read sequences without end-to-end alignments
        overlap. In general, readlets have size max_readlet_size and are
        spaced readlet_interval bases apart. Capping readlets cover the ends of
        reads and have minimum size min_readlet_size. capping_multiplier is 
        used to obtain the lengths of other capping readlets according to the 
        formula min_readlet_size*capping_multiplier^n for n an integer >= 0.

        seq: sequence to readletize
        sample_indexes: dictionary mapping indexes of samples in which seq
            appears to their counts in that sample
        reversed_complement_sample_indexes: dictionary mapping indexes of
            samples in which reversed complement of seq appears to their counts
            in that sample
        seq_id: unique identifier to assign to a given sequence
        cap_sizes: list of sizes of capping readlets
        min_readlet_size: "capping" readlets (that is, readlets that terminate
            at a given end of the read) are never smaller than this value.
            Ignored if readletize=False.
        max_readlet_size: size of every noncapping readlet. Ignored if 
            readletize=False.
        readlet_interval: number of bases separating successive readlets along
            the read
        verbose: True if reads/readlets should occasionally be written 
            to stderr.
        report_multiplier: if verbose is True, the line number of a read or its 
            first readlet written to stderr increases exponentially with base
            report_multiplier.
        no_polyA: kill readlets that are all As

        No return value.
    """
    global _output_line_count
    if len(seq) < min_readlet_size:
        return
    '''Construct a readlet identifier as follows: read sequence ID + 
    ('-' if readlet sequence is reverse-complemented; else '+') + '\x1e' +
    displacement of readlet's 5' end from read's 5' end + '\x1e' +
    displacement of readlet's 3' end from read's 3' end (+, for EXACTLY one
    readlet of a read sequence, '\x1e' + read sequence + '\x1e' + number of
    instances of read sequence + '\x1e' + number of instances of read
    sequence's reversed complement + '\x1e' + (an '\x1f'-separated list A of
    unique sample labels with read sequences that match the original read
    sequence) + '\x1e' + (an '\x1f'-separated list  of unique sample labels B
    with read sequences that match the reversed complement of the original read
    sequence)) + '\x1e' + (an '\x1f'-separated list of the number of instances
    of the read sequence for each respective sample in list A) + '\x1e' + (an
    '\x1f'-separated list of the number of instances of the read
    sequence's reversed complement for each respective sample in list B). Here,
    a read sequence ID takes the form X:Y, where X is the
    "mapred_task_partition" environment variable -- a unique index for a task
    within a job -- and Y is the index of the read sequence relative to the
    beginning of the input stream.'''
    to_write = []
    seq_size = len(seq)
    # Add capping readlets
    for cap_size in cap_sizes:
        readlet_seq = seq[:cap_size]
        reversed_complement_readlet_seq = readlet_seq[::-1].translate(
                        _reversed_complement_translation_table
                    )
        if readlet_seq < reversed_complement_readlet_seq:
            if not no_polyA or set(readlet_seq) != _polyA:
                to_write.append('%s\t%s\x1e%d\x1e%d' % (readlet_seq, seq_id 
                                                        + '+', 0,
                                                        seq_size - cap_size))
        else:
            if not no_polyA or set(reversed_complement_readlet_seq) != _polyA:
                to_write.append('%s\t%s\x1e%d\x1e%d' \
                                    % (reversed_complement_readlet_seq,
                                        seq_id + '-', 0, seq_size - cap_size))
        readlet_seq = seq[-cap_size:]
        reversed_complement_readlet_seq = readlet_seq[::-1].translate(
                        _reversed_complement_translation_table
                    )
        if readlet_seq < reversed_complement_readlet_seq:
            if not no_polyA or set(readlet_seq) != _polyA:
                to_write.append('%s\t%s\x1e%d\x1e%d' % (readlet_seq,
                                                    seq_id + '+',
                                                    seq_size - cap_size, 0))
        else:
            if not no_polyA or set(reversed_complement_readlet_seq) != _polyA:
                to_write.append('%s\t%s\x1e%d\x1e%d' \
                                    % (reversed_complement_readlet_seq,
                                        seq_id + '-', seq_size - cap_size, 0))
    # Add noncapping readlets
    for j in xrange(readlet_interval, seq_size - max_readlet_size,
                        readlet_interval):
        readlet_seq = seq[j:j+max_readlet_size]
        reversed_complement_readlet_seq = readlet_seq[::-1].translate(
                        _reversed_complement_translation_table
                    )
        if readlet_seq < reversed_complement_readlet_seq:
            if not no_polyA or set(readlet_seq) != _polyA:
                to_write.append('%s\t%s\x1e%d\x1e%d' % (readlet_seq,
                                                    seq_id + '+', j, 
                                                    seq_size - j
                                                        - max_readlet_size))
        else:
            if not no_polyA or set(reversed_complement_readlet_seq) != _polyA:
                to_write.append('%s\t%s\x1e%d\x1e%d' \
                                    % (reversed_complement_readlet_seq,
                                        seq_id + '-', j, 
                                        seq_size - j - max_readlet_size))
    # Add additional info to first readlet in to_write
    try:
        to_write[0] = '\x1e'.join([to_write[0], seq,
                        '\x1f'.join(sample_indexes.keys()),
                        '\x1f'.join(
                            reversed_complement_sample_indexes.keys()
                        ),
                        '\x1f'.join(
                            [str(count) for count in sample_indexes.values()]
                        ),
                        '\x1f'.join(
                            [str(count) for count
                                in reversed_complement_sample_indexes.values()]
                        )])
    except IndexError:
        # No readlets because polyAs were present!
        return
    if verbose:
        print >>sys.stderr, 'First readlet from read %d: %s' \
            % (i + 1, to_write[0])
        next_report_line = int((next_report_line + 1)
            * report_multiplier + 1) - 1
    for readlet in to_write:
        print >>output_stream, 'readletized\t%s' % readlet
    _output_line_count += len(to_write)

def qname_and_mate(qname):
    """ Removes mate sequence from qname.

        qname: original qname

        Return value: tuple (qname, mate sequence, which could be empty)
    """
    split_name = qname.split('\x1d')
    mate = split_name[1].rpartition(':')
    split_name[1] = mate[0] + ':'
    return ('\x1d'.join(split_name), mate[2])

def handle_bowtie_output(input_stream, reference_index, manifest_object,
        group_reads_object, task_partition, cap_sizes, k_value=1,
        align_stream=None, other_stream=None, output_stream=sys.stdout,
        exon_differentials=True, exon_intervals=False, verbose=False,
        bin_size=10000, report_multiplier=1.2, search_filter=8,
        min_readlet_size=8, max_readlet_size=25,
        readlet_interval=5, drop_deletions=False, output_bam_by_chr=False,
        tie_margin=0, no_realign=False, no_polyA=False):
    """ Prints end-to-end alignments and selects reads to be realigned.

        input_stream: where to retrieve Bowtie's SAM output, typically a
            Bowtie process's stdout.
        reference_index: object of class bowtie_index.BowtieIndexReference
            that permits access to reference
        manifest_object: object of class LabelsAndIndices that maps indices
            to labels and back; used to shorten intermediate output.
        group_reads_object: object of class IndexGroup for grapping
            index group for a given read sequence
        task_partition: unique identifier for task; derived from environment
            variable and used for form sequence identifier for unique read
            sequence
        cap_sizes: list of sizes of capping readlets
        k_value: argument of Bowtie 2's -k parameter
        align_stream: where to write reads that are to be aligned on second
            pass
        other_stream: where reads with the same sequence as those from the
            input stream are stored; if the primary alignment of a read from
            the input stream is not a tie and the edit distance (NM:i) is 0,
            these reads are assigned the same alignments. If None, second-pass
            alignment is being performed because there were ties or a nonzero
            edit distance during first-pass alignment.
        output_stream: where to emit exon and junction tuples; typically,
            this is sys.stdout.
        exon_differentials: True iff EC differentials are to be emitted.
        exon_intervals: True iff EC intervals are to be emitted.
        verbose: True if alignments should occasionally be written 
            to stderr.
        bin_size: genome is partitioned in units of bin_size for later load
            balancing.
        report_multiplier: if verbose is True, the line number of an
            alignment written to stderr increases exponentially with base
            report_multiplier.
        search_filter: determine how large a soft clip on one side of a read
            is necessary to pass it on to junction search pipeline
        min_readlet_size: "capping" readlets (that is, readlets that terminate
            at a given end of the read) are never smaller than this value.
            Ignored if readletize=False.
        max_readlet_size: size of every noncapping readlet. Ignored if 
            readletize=False.
        readlet_interval: number of bases separating successive readlets along
            the read
        drop_deletions: True iff deletions should be dropped from coverage
            vector
        output_bam_by_chr: True iff final BAMs are to be output by chromosome
        tie_margin: allowed score difference per 100 bases among ties in
            max score. For example, 150 and 144 are tied alignment scores
            for a 100-bp read when --tie-margin is 6
        no_realign: True iff job flow does not need more than readlets: this
            usually means only a transcript index is being constructed
        no_polyA: kill readlets that are all As

        No return value.
    """
    global _output_line_count
    next_report_line = 0
    i = 0
    # Shortcut if no realignment is to be done in job flow
    if other_stream and no_realign:
        # First-pass alignment
        if verbose:
            print >>sys.stderr, ('Processing alignments; suppressing some '
                                 'output.')
        readletized_index = 0
        other_xstream = xstream(other_stream, 1)
        for (qname,), xpartition in xstream(input_stream, 1):
            is_reverse, _, qname = qname.partition('\x1d')
            qname, _ = qname_and_mate(qname)
            is_reverse = int(is_reverse)
            _, other_xpartition = other_xstream.next()
            sample_index = manifest_object.label_to_index[
                        qname.rpartition('\x1d')[2]
                    ]
            # Handle primary alignment
            rest_of_line = xpartition.next()
            i += 1
            flag = int(rest_of_line[0])
            if is_reverse:
                true_flag = flag ^ 16
            else:
                true_flag = flag
            cigar = rest_of_line[4]
            rname = rest_of_line[1]
            pos = int(rest_of_line[2])
            seq, qual = rest_of_line[8], rest_of_line[9]
            if verbose and next_report_line == i:
                print >>sys.stderr, \
                    'SAM output record %d: rdname="%s", flag=%d' \
                    % (i, qname, true_flag)
                next_report_line = max(int(next_report_line
                    * report_multiplier), next_report_line + 1)
            if flag & 4:
                continue
            multiread = [(qname,) + rest_of_line] + \
                [(qname,) + next_line for next_line in xpartition]
            scores = [
                    field[5:] for field in multiread[0]
                    if field[:5] == 'XS:i:'
                ] + [
                    field[5:] for field in multiread[0]
                    if field[:5] == 'AS:i:'
                ]
            # Note: needs to be changed to accommodate --tie-margin
            # extra keywords for grepping: tie margin tie_margin break
            tie_present = (len(scores) == 2) and (scores[0] == scores[1])
            clip_present = 'S' in cigar
            exact_match = 'NM:i:0' in multiread[0]
            alignment_reversed = [
                                    int(alignment[1]) & 16
                                    for alignment in multiread
                                ]
            if k_value == 1 and exact_match and not clip_present:
                continue
            reversed_complement_seq = seq[::-1].translate(
                _reversed_complement_translation_table
            )
            if seq < reversed_complement_seq:
                seq_to_print = seq
                qual_to_print = qual
            else:
                seq_to_print = reversed_complement_seq
                qual_to_print = qual[::-1]
            if not (exact_match and not clip_present):
                # Print output
                split_cigar = re.split(r'([MINDS])', cigar)[:-1]
                try:
                    if ((split_cigar[1] == 'S'
                            and int(split_cigar[0]) >= search_filter) or
                        (split_cigar[-1] == 'S'
                            and int(split_cigar[-2]) >= search_filter)):
                        search_for_junctions = True
                    else:
                        search_for_junctions = False
                except IndexError:
                    search_for_junctions = False
                if search_for_junctions:
                    sample_indexes, reversed_complement_sample_indexes = (
                            defaultdict(int), defaultdict(int)
                        )
                    try:
                        for current_is_reverse, current_qname, current_qual \
                            in other_xpartition:
                            if current_is_reverse == '1':
                                reversed_complement_sample_indexes[
                                        manifest_object.label_to_index[
                                            current_qname.rpartition('\x1d')[2]
                                        ]
                                    ] += 1
                            else:
                                sample_indexes[
                                        manifest_object.label_to_index[
                                            current_qname.rpartition('\x1d')[2]
                                        ]
                                    ] += 1
                    except ValueError:
                        # No other reads
                        pass
                    if is_reverse:
                        reversed_complement_sample_indexes[sample_index] += 1
                    else:
                        sample_indexes[sample_index] += 1
                    print_readletized_output(
                            seq=seq_to_print,
                            sample_indexes=sample_indexes,
                            reversed_complement_sample_indexes=\
                                reversed_complement_sample_indexes,
                            seq_id=(':'.join([task_partition,
                                                str(readletized_index)])),
                            cap_sizes=cap_sizes,
                            output_stream=output_stream,
                            min_readlet_size=min_readlet_size,
                            max_readlet_size=max_readlet_size,
                            readlet_interval=readlet_interval,
                            verbose=(verbose and next_report_line == i),
                            no_polyA=no_polyA)
                    readletized_index += 1
        return
    alignment_printer = AlignmentPrinter(
            manifest_object,
            reference_index,
            bin_size=bin_size,
            output_stream=output_stream,
            exon_ivals=exon_intervals,
            exon_diffs=exon_differentials,
            drop_deletions=drop_deletions,
            output_bam_by_chr=output_bam_by_chr,
            tie_margin=tie_margin
        )
    if other_stream:
        # First-pass alignment
        if verbose:
            print >>sys.stderr, 'Processing first-pass alignments.'
        readletized_index = 0
        other_xstream = xstream(other_stream, 1)
        for (qname,), xpartition in xstream(input_stream, 1):
            is_reverse, _, qname = qname.partition('\x1d')
            qname, mate = qname_and_mate(qname)
            is_reverse = int(is_reverse)
            _, other_xpartition = other_xstream.next()
            sample_index = manifest_object.label_to_index[
                        qname.rpartition('\x1d')[2]
                    ]
            # Handle primary alignment
            rest_of_line = xpartition.next()
            i += 1
            flag = int(rest_of_line[0])
            if is_reverse:
                true_flag = flag ^ 16
            else:
                true_flag = flag
            cigar = rest_of_line[4]
            rname = rest_of_line[1]
            pos = int(rest_of_line[2])
            seq, qual = rest_of_line[8], rest_of_line[9]
            if verbose and next_report_line == i:
                print >>sys.stderr, \
                    'SAM output record %d: rdname="%s", flag=%d' \
                    % (i, qname, true_flag)
                next_report_line = max(int(next_report_line
                    * report_multiplier), next_report_line + 1)
            if flag & 4:
                index_partition = group_reads_object.index_group(seq)
                print >>output_stream, 'unique\t%s' % seq
                print >>output_stream, 'unmapped\t%s\t%s\t%d\t%s\t%s' % (
                                                            index_partition,
                                                            seq,
                                                            is_reverse + 1,
                                                            qname,
                                                            qual
                                                        )
                try:
                    for is_reverse, qname, qual in other_xpartition:
                        qname, _ = qname_and_mate(qname)
                        print >>output_stream, (
                                        'unmapped\t%s\t%s\t%d\t%s\t%s' % (
                                                            index_partition,
                                                            seq,
                                                            int(is_reverse)
                                                                + 1,
                                                            qname,
                                                            qual
                                                        )
                                    )
                except ValueError:
                    pass
                continue
            multiread = [(qname,) + rest_of_line] + \
                [(qname,) + next_line for next_line in xpartition]
            scores = [
                    field[5:] for field in multiread[0]
                    if field[:5] == 'XS:i:'
                ] + [
                    field[5:] for field in multiread[0]
                    if field[:5] == 'AS:i:'
                ]
            # Note: needs to be changed to accommodate --tie-margin
            # extra keywords for grepping: tie margin tie_margin break
            tie_present = (len(scores) == 2) and (scores[0] == scores[1])
            clip_present = 'S' in cigar
            exact_match = 'NM:i:0' in multiread[0]
            alignment_reversed = [
                                    int(alignment[1]) & 16
                                    for alignment in multiread
                                ]
            if k_value == 1 and exact_match and not clip_present:
                NH_field = 'NH:i:%d' % len(multiread)
                reversed_flags = None
                if not tie_present:
                    # All alignments for all read sequences can be written
                    try:
                        for current_is_reverse, current_qname, current_qual \
                            in other_xpartition:
                            if current_is_reverse == '1':
                                if not reversed_flags:
                                    reversed_flags = [
                                            str(int(alignment[1]) ^ 16)
                                            for alignment in multiread
                                        ]
                                _output_line_count += \
                                    alignment_printer.print_alignment_data(
                                            ([(current_qname,
                                                reversed_flags[j])
                                               + alignment[2:4]
                                               + ('255',)
                                               + alignment[5:10]
                                               + (current_qual[::-1]
                                                    if alignment_reversed[j]
                                                    else current_qual,)
                                               + tuple([
                                                    field for field in
                                                    alignment[11:] if 
                                                    field[:5] != 'XS:i:'
                                                ])
                                               + (NH_field,)
                                               for j, alignment
                                               in enumerate(multiread)],)
                                        )
                            else:
                                _output_line_count += \
                                    alignment_printer.print_alignment_data(
                                            ([(current_qname,)
                                                + alignment[1:4]
                                                + ('255',)
                                                + alignment[5:10]
                                                + (current_qual[::-1]
                                                  if alignment_reversed[j]
                                                  else current_qual,)
                                                + tuple([
                                                    field for field in
                                                    alignment[11:] if 
                                                    field[:5] != 'XS:i:'
                                                ])
                                               + (NH_field,)
                                               for j, alignment
                                               in enumerate(multiread)],)
                                        )
                    except ValueError:
                        pass
                if tie_present and mate:
                    # Write read _pair_ to align_stream
                    decoded = decode_sequence(mate)
                    if is_reverse:
                        if flag & 16:
                            print >>align_stream, '\t'.join([
                                    qname, seq, qual,
                                        decoded, len(decoded)*'I'
                                ])
                        else:
                            print >>align_stream, '\t'.join([
                                    qname, seq[::-1].translate(
                                        _reversed_complement_translation_table
                                    ), qual[::-1], decoded, len(decoded)*'I'
                                ])
                    else:
                        if flag & 16:
                            print >>align_stream, '\t'.join([
                                    qname, seq[::-1].translate(
                                        _reversed_complement_translation_table
                                    ), qual[::-1], decoded, len(decoded)*'I'
                                ])
                        else:
                            print >>align_stream, '\t'.join([
                                    qname, seq, qual, decoded, len(decoded)*'I'
                                ])
                else:
                    # Final alignment can be written if no tie+mate
                    if is_reverse:
                        _output_line_count += \
                            alignment_printer.print_alignment_data(
                                    ([(alignment[0],
                                        str(int(alignment[1]) ^ 16))
                                       + alignment[2:] + (NH_field,)
                                       for alignment in multiread],)
                                )
                    else:
                        _output_line_count += \
                            alignment_printer.print_alignment_data(
                                    ([alignment + (NH_field,)
                                       for alignment in multiread],)
                                )
            else:
                print >>output_stream, 'unique\t%s' % seq
            reversed_complement_seq = seq[::-1].translate(
                _reversed_complement_translation_table
            )
            if seq < reversed_complement_seq:
                seq_to_print = seq
                qual_to_print = qual
            else:
                seq_to_print = reversed_complement_seq
                qual_to_print = qual[::-1]
            if flag & 16:
                seq_to_realign = reversed_complement_seq
                reversed_complement_seq_to_realign = seq
                qual_to_realign = qual[::-1]
                reversed_qual_to_realign = qual
            else:
                seq_to_realign = seq
                reversed_complement_seq_to_realign = reversed_complement_seq
                qual_to_realign = qual
                reversed_qual_to_realign = qual[::-1]
            if not (exact_match and not clip_present):
                # Prep for second-pass Bowtie 2 and print output
                split_cigar = re.split(r'([MINDS])', cigar)[:-1]
                try:
                    if ((split_cigar[1] == 'S'
                            and int(split_cigar[0]) >= search_filter) or
                        (split_cigar[-1] == 'S'
                            and int(split_cigar[-2]) >= search_filter)):
                        search_for_junctions = True
                    else:
                        search_for_junctions = False
                except IndexError:
                    search_for_junctions = False
                if search_for_junctions:
                    sample_indexes, reversed_complement_sample_indexes = (
                            defaultdict(int), defaultdict(int)
                        )
                    try:
                        for current_is_reverse, current_qname, current_qual \
                            in other_xpartition:
                            current_qname, current_mate = qname_and_mate(
                                    current_qname
                                )
                            if tie_present and current_mate:
                                decoded = decode_sequence(current_mate)
                                if current_is_reverse == '1':
                                    print >>align_stream, '\t'.join([
                                            current_qname,
                                            reversed_complement_seq_to_realign,
                                            current_qual[::-1],
                                            decoded, 'I'*len(decoded)
                                    ])
                                else:
                                    print >>align_stream, '\t'.join([
                                            current_qname,
                                            seq_to_realign, current_qual,
                                            decoded, 'I'*len(decoded)
                                    ])
                            else:
                                if current_is_reverse == '1':
                                    print >>align_stream, '\t'.join([
                                            current_qname,
                                            reversed_complement_seq_to_realign,
                                            current_qual[::-1]
                                        ])
                                else:
                                    print >>align_stream, '\t'.join([
                                            current_qname, seq_to_realign,
                                            current_qual
                                        ])
                            if current_is_reverse == '1':
                                reversed_complement_sample_indexes[
                                        manifest_object.label_to_index[
                                            current_qname.rpartition('\x1d')[2]
                                        ]
                                    ] += 1
                            else:
                                sample_indexes[
                                        manifest_object.label_to_index[
                                            current_qname.rpartition('\x1d')[2]
                                        ]
                                    ] += 1
                    except ValueError:
                        # No other reads
                        pass
                    if is_reverse:
                        reversed_complement_sample_indexes[sample_index] += 1
                    else:
                        sample_indexes[sample_index] += 1
                    print_readletized_output(
                            seq=seq_to_print,
                            sample_indexes=sample_indexes,
                            reversed_complement_sample_indexes=\
                                reversed_complement_sample_indexes,
                            seq_id=(':'.join([task_partition,
                                                str(readletized_index)])),
                            cap_sizes=cap_sizes,
                            output_stream=output_stream,
                            min_readlet_size=min_readlet_size,
                            max_readlet_size=max_readlet_size,
                            readlet_interval=readlet_interval,
                            verbose=(verbose and next_report_line == i),
                            no_polyA=no_polyA)
                    readletized_index += 1
                elif tie_present:
                    try:
                        for current_is_reverse, current_qname, current_qual \
                            in other_xpartition:
                            current_qname, current_mate = qname_and_mate(
                                    current_qname
                                )
                            if current_mate:
                                decoded = decode_sequence(current_mate)
                                if current_is_reverse == '1':
                                    print >>align_stream, '\t'.join([
                                            current_qname,
                                            reversed_complement_seq_to_realign,
                                            current_qual[::-1],
                                            decoded, 'I'*len(decoded)
                                    ])
                                else:
                                    print >>align_stream, '\t'.join([
                                            current_qname,
                                            seq_to_realign, current_qual,
                                            decoded, 'I'*len(decoded)
                                    ])
                            else:
                                if current_is_reverse == '1':
                                    print >>align_stream, '\t'.join([
                                            current_qname,
                                            reversed_complement_seq_to_realign,
                                            current_qual[::-1]
                                        ])
                                else:
                                    print >>align_stream, '\t'.join([
                                            current_qname, seq_to_realign,
                                            current_qual
                                        ])
                    except ValueError:
                        # No other reads
                        pass
                else:
                    try:
                        for current_is_reverse, current_qname, current_qual \
                            in other_xpartition:
                            current_qname, _ = qname_and_mate(
                                    current_qname
                                )
                            if current_is_reverse == '1':
                                print >>align_stream, '\t'.join([
                                        current_qname,
                                        reversed_complement_seq_to_realign,
                                        current_qual[::-1]
                                    ])
                            else:
                                print >>align_stream, '\t'.join([
                                        current_qname, seq_to_realign,
                                        current_qual
                                    ])
                    except ValueError:
                        # No other reads
                        pass
            elif tie_present:
                # Prep for second-pass Bowtie 2 with (possibly) pairs
                try:
                    for current_is_reverse, current_qname, current_qual \
                        in other_xpartition:
                        current_qname, current_mate = qname_and_mate(
                                current_qname
                            )
                        if current_mate:
                            decoded = decode_sequence(current_mate)
                            if current_is_reverse == '1':
                                print >>align_stream, '\t'.join([
                                        current_qname,
                                        reversed_complement_seq_to_realign,
                                        current_qual[::-1],
                                        decoded, 'I'*len(decoded)
                                ])
                            else:
                                print >>align_stream, '\t'.join([
                                        current_qname,
                                        seq_to_realign, current_qual,
                                        decoded, 'I'*len(decoded)
                                ])
                        else:
                            if current_is_reverse == '1':
                                print >>align_stream, '\t'.join([
                                        current_qname,
                                        reversed_complement_seq_to_realign,
                                        current_qual[::-1]
                                    ])
                            else:
                                print >>align_stream, '\t'.join([
                                        current_qname, seq_to_realign,
                                        current_qual
                                    ])
                except ValueError:
                    # No other reads
                    pass
            elif k_value != 1:
                # Only prepare for second-pass Bowtie 2
                try:
                    for current_is_reverse, current_qname, current_qual \
                        in other_xpartition:
                        current_qname, _ = qname_and_mate(current_qname)
                        if current_is_reverse == '1':
                            print >>align_stream, '\t'.join([
                                    current_qname,
                                    reversed_complement_seq_to_realign,
                                    current_qual[::-1]
                                ])
                        else:
                            print >>align_stream, '\t'.join([
                                    current_qname, seq_to_realign,
                                    current_qual
                                ])
                except ValueError:
                    # No other reads
                    pass
            if k_value != 1 or not (exact_match and not clip_present):
                '''Write "postponed" SAM/unmapped lines if there's no tie+mate;
                otherwise, try to align again.'''
                if tie_present and mate:
                    decoded = decode_sequence(mate)
                    if is_reverse:
                        print >>align_stream, '\t'.join([
                                        qname,
                                        reversed_complement_seq_to_realign,
                                        reversed_qual_to_realign,
                                        decoded, 'I'*len(decoded)
                                    ])
                    else:
                        print >>align_stream, '\t'.join([
                                        qname, seq_to_realign, qual_to_realign,
                                        decoded, 'I'*len(decoded)
                                    ])
                else:
                    if is_reverse:
                        for alignment in multiread:
                            print >>output_stream, \
                                '\t'.join(('postponed_sam', alignment[0])
                                    + (str(int(alignment[1]) ^ 16),)
                                    + alignment[2:])
                            _output_line_count += 1
                    else:
                        for alignment in multiread:
                            print >>output_stream, \
                                '\t'.join(('postponed_sam',) + alignment)
                            _output_line_count += 1
                    print >>output_stream, 'unmapped\t%s\t%s\t%d\t%s\t%s' % (
                                group_reads_object.index_group(seq_to_print),
                                seq_to_print,
                                is_reverse + 1,
                                qname,
                                qual_to_print
                            )
    else:
        # Second-pass alignment
        if verbose:
            print >>sys.stderr, 'Processing second-pass alignments.'
        for (qname,), xpartition in xstream(input_stream, 1):
            # Add to multiread if and only if not second mate from pair
            multiread = [(qname,) + next_line for next_line in xpartition
                            if not (int(next_line[0]) & 128)]
            # Handle primary alignment
            rest_of_line = multiread[0][1:]
            i += 1
            flag = int(rest_of_line[0])
            cigar = rest_of_line[4]
            rname = rest_of_line[1]
            pos = int(rest_of_line[2])
            seq, qual = rest_of_line[8], rest_of_line[9]
            reversed_complement_seq = seq[::-1].translate(
                _reversed_complement_translation_table
            )
            if seq < reversed_complement_seq:
                seq_to_print = seq
                qual_to_print = qual
                if flag & 16:
                    is_reverse = 1
                else:
                    is_reverse = 0
            else:
                seq_to_print = reversed_complement_seq
                qual_to_print = qual[::-1]
                if flag & 16:
                    is_reverse = 0
                else:
                    is_reverse = 1
            if verbose and next_report_line == i:
                print >>sys.stderr, \
                    'SAM output record %d: rdname="%s", flag=%d' \
                    % (i, qname, flag)
                next_report_line = max(int(next_report_line
                    * report_multiplier), next_report_line + 1)
            if flag & 4:
                print >>output_stream, 'unmapped\t%s\t%s\t%d\t%s\t%s' % (
                                group_reads_object.index_group(seq_to_print),
                                seq_to_print,
                                is_reverse + 1,
                                qname,
                                qual_to_print
                            )
                continue
            clip_present = 'S' in cigar
            exact_match = 'NM:i:0' in multiread[0]
            if k_value == 1 and exact_match and not clip_present:
                # All alignments for all read sequences can be written
                NH_field = 'NH:i:%d' % len(multiread)
                _output_line_count += \
                    alignment_printer.print_alignment_data(
                            ([alignment + (NH_field,)
                               for alignment in multiread],)
                        )
            else:
                # Write "postponed" SAM/unmapped lines
                for alignment in multiread:
                    print >>output_stream, \
                        '\t'.join(('postponed_sam',) + alignment)
                    _output_line_count += 1
                reversed_complement_seq = seq[::-1].translate(
                    _reversed_complement_translation_table
                )
                print >>output_stream, 'unmapped\t%s\t%s\t%d\t%s\t%s' % (
                                group_reads_object.index_group(seq_to_print),
                                seq_to_print,
                                is_reverse + 1,
                                qname,
                                qual_to_print
                            )
    output_stream.flush()

def go(task_partition='0', other_reads=None, second_pass_reads=None,
        min_readlet_size=15, drop_deletions=False, max_readlet_size=25,
        readlet_interval=4, capping_multiplier=1.5, output_stream=sys.stdout,
        input_stream=sys.stdin, verbose=False, report_multiplier=1.2,
        k_value=1, bowtie_index_base='genome', bin_size=10000,
        manifest_file='manifest', exon_differentials=True,
        exon_intervals=False, gzip_level=3, search_filter=9,
        index_count=1, output_bam_by_chr=False, tie_margin=0,
        no_realign=False, no_polyA=False):
    """ Emits output specified in align_reads.py by processing Bowtie 2 output.

        This script containing this function is invoked twice to process each
        of two passes of Bowtie 2. On first pass, every unique read sequence
        input to the step is aligned exactly once. The quality string for this
        sequence has the highest mean quality from among the set of reads that
        share the same sequence. (See align_reads.py for details on how the
        quality string is obtained.) If exactly one perfect match of the read
        sequence is found in the reference, that match is assigned to all reads
        with the same sequence. If a first-pass alignment is soft-clipped, the
        corresponding read sequence is readletized, and other reads with the
        same sequence are not realigned on second pass. Otherwise, read
        sequences that share the same sequence are realigned on second-pass
        Bowtie 2. All imperfect alignments are saved for realignment to
        transcript fragments in later steps.

        input_stream: where to retrieve Bowtie 2 output
        output_stream: where to emit output specified in align_reads.py
        task_partition: MapReduce task partition; used to assign an identifier
            to each unique read sequence
        other_reads: where to find reads that are not to be aligned on
            first pass but may be aligned on second pass; None if second
            invocation of script
        second_pass_reads: file to which reads that are to be aligned on 
            second pass should be written; None if second invocation of script
        bowtie_index_base: the basename of the Bowtie1 index files associated
            with the reference.
        k_value: argument of -k from Bowtie 2 arguments
        manifest_file: filename of manifest
        bin_size: genome is partitioned in units of bin_size for later load
            balancing
        verbose: True if alignments should occasionally be written to stderr.
        exon_differentials: True iff EC differentials are to be emitted.
        exon_intervals: True iff EC intervals are to be emitted.
        report_multiplier: if verbose is True, the line number of an alignment
            or read written to stderr increases exponentially with base
            report_multiplier.
        search_filter: determines how large a soft clip on one side of
            a read is necessary to pass it on to junction search pipeline
        min_readlet_size: "capping" readlets (that is, readlets that terminate
            at a given end of the read) are never smaller than this value
        max_readlet_size: size of every noncapping readlet
        readlet_interval: number of bases separating successive readlets along
            the read
        capping_multiplier: successive capping readlets on a given end of a
            read are increased in size exponentially with base
            capping_multiplier
        drop_deletions: True iff deletions should be dropped from coverage
            vector
        gzip_level: gzip level to use to compress temporary files
        index_count: number of transcriptome Bowtie 2 indexes to which to
            assign unmapped reads for later realignment
        output_bam_by_chr: True iff final BAM output should be by chr
        tie_margin: allowed score difference per 100 bases among ties in
            max score. For example, 150 and 144 are tied alignment scores
            for a 100-bp read when --tie-margin is 6
        no_realign: True iff job flow does not need more than readlets: this
            usually means only a transcript index is being constructed
        no_polyA: kill readlets that are all As
    """
    reference_index = bowtie_index.BowtieIndexReference(bowtie_index_base)
    manifest_object = manifest.LabelsAndIndices(manifest_file)
    group_reads_object = group_reads.IndexGroup(index_count)
    if other_reads is not None:
        # First-pass Bowtie 2
        assert second_pass_reads is not None
        # Build list of capping readlet sizes
        cap_sizes = []
        cap_size = min_readlet_size
        while cap_size <= max_readlet_size:
            cap_sizes.append(cap_size)
            cap_size = int(cap_size*capping_multiplier)
            if cap_size == cap_sizes[-1]:
              cap_size += 1
        if cap_size != max_readlet_size:
            # Always have a start or end read of length max_readlet_size
            cap_sizes.append(max_readlet_size)
        with xopen(None, other_reads) as other_stream, \
            xopen(True, second_pass_reads, 'w') as align_stream:
            handle_bowtie_output(
                    input_stream,
                    reference_index,
                    manifest_object,
                    group_reads_object,
                    task_partition,
                    cap_sizes,
                    k_value=k_value,
                    align_stream=align_stream,
                    other_stream=other_stream,
                    exon_differentials=exon_differentials, 
                    exon_intervals=exon_intervals, 
                    bin_size=bin_size,
                    verbose=verbose, 
                    output_stream=output_stream,
                    report_multiplier=report_multiplier,
                    search_filter=search_filter,
                    min_readlet_size=min_readlet_size,
                    max_readlet_size=max_readlet_size,
                    readlet_interval=readlet_interval,
                    drop_deletions=drop_deletions,
                    output_bam_by_chr=output_bam_by_chr,
                    tie_margin=tie_margin,
                    no_realign=no_realign,
                    no_polyA=no_polyA
                )
        print >>sys.stderr, (
            'align_reads_delegate.py reports %d output lines on first pass.'
            % _output_line_count
        )
    else:
        handle_bowtie_output(
                input_stream,
                reference_index,
                manifest_object,
                group_reads_object,
                task_partition,
                [],
                k_value=k_value,
                exon_differentials=exon_differentials, 
                exon_intervals=exon_intervals, 
                bin_size=bin_size,
                verbose=verbose, 
                output_stream=output_stream,
                report_multiplier=report_multiplier,
                search_filter=search_filter,
                drop_deletions=drop_deletions,
                output_bam_by_chr=output_bam_by_chr,
                tie_margin=tie_margin
            )
        print >>sys.stderr, (
            'align_reads_delegate.py reports %d output lines on second pass.'
            % _output_line_count
        )
    output_stream.flush()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--verbose', action='store_const', const=True,
        default=False,
        help='Print out extra debugging statements')
    parser.add_argument('--task-partition', type=str, required=False,
        default=None,
        help='Task partition ID; used to label unique read sequences')
    parser.add_argument('--other-reads', type=str, required=False,
        default=None,
        help=('Path to file containing reads that are not to be aligned on '
              'first pass; included only on first invocation of script'))
    parser.add_argument('--second-pass-reads', type=str, required=False,
        default=None,
        help=('Path to file to write containing reads that are to be aligned '
              'on second pass; included only on first invocation of script'))
    parser.add_argument('--drop-deletions', action='store_const',
        const=True,
        default=False, 
        help='Drop deletions from coverage vectors')
    parser.add_argument('--output-bam-by-chr', action='store_const',
        const=True,
        default=False, 
        help='Final BAMs will be output by chromosome')
    parser.add_argument('--exon-differentials', action='store_const',
        const=True,
        default=True, 
        help='Print exon differentials (+1s and -1s)')
    parser.add_argument('--exon-intervals', action='store_const',
        const=True,
        default=False, 
        help='Print exon intervals')
    parser.add_argument('--search-filter', type=int, required=False,
        default=1,
        help=('Minimum size of soft-clipped end that dispatches a read for '
              'junction search'))
    parser.add_argument('--k-value', type=int, required=False,
        default=1,
        help='Argument of -k parameter from Bowtie 2')
    parser.add_argument('--min-readlet-size', type=int, required=False,
        default=15, 
        help='Capping readlets (that is, readlets that terminate '
             'at a given end of the read) are never smaller than this value')
    parser.add_argument('--max-readlet-size', type=int, required=False,
        default=25, 
        help='Size of every noncapping readlet')
    parser.add_argument('--readlet-interval', type=int, required=False,
        default=12, 
        help='Number of bases separating successive noncapping readlets along '
             'the read')
    parser.add_argument('--capping-multiplier', type=float, required=False,
        default=1.5, 
        help='Successive capping readlets on a given end of a read are '
             'increased in size exponentially with this base')
    parser.add_argument('--report-multiplier', type=float, required=False,
        default=1.2,
        help='When --verbose is also invoked, the only lines of lengthy '
             'intermediate output written to stderr have line number that '
             'increases exponentially with this base')
    parser.add_argument('--test', action='store_const', const=True,
        default=False,
        help='Run unit tests; DOES NOT NEED INPUT FROM STDIN')
    parser.add_argument('--gzip-level', type=int, required=False,
        default=3,
        help='Level of gzip compression to use for temporary files')
    parser.add_argument('--no-realign', action='store_const',
        const=True,
        default=False, 
        help=('Suppresses some output and does not perform second-pass '
              'alignment if unnecessary for job flow'))
    parser.add_argument('--no-polyA', action='store_const',
        const=True,
        default=False, 
        help='Disallows any readlet that is a string of A nucleotides')

    # Add command-line arguments for dependencies
    partition.add_args(parser)
    bowtie.add_args(parser)
    group_reads.add_args(parser)
    manifest.add_args(parser)
    from alignment_handlers import add_args as alignment_handlers_add_args
    alignment_handlers_add_args(parser)

    args = parser.parse_args()

if __name__ == '__main__' and not args.test:
    go(task_partition=args.task_partition,
        other_reads=args.other_reads,
        second_pass_reads=args.second_pass_reads,
        min_readlet_size=args.min_readlet_size,
        drop_deletions=args.drop_deletions,
        max_readlet_size=args.max_readlet_size,
        readlet_interval=args.readlet_interval,
        capping_multiplier=args.capping_multiplier,
        output_stream=sys.stdout,
        input_stream=sys.stdin,
        verbose=args.verbose,
        report_multiplier=args.report_multiplier,
        k_value=args.k_value,
        bowtie_index_base=args.bowtie_idx,
        bin_size=args.partition_length,
        manifest_file=args.manifest,
        exon_differentials=args.exon_differentials,
        exon_intervals=args.exon_intervals,
        search_filter=args.search_filter,
        gzip_level=args.gzip_level,
        index_count=args.index_count,
        output_bam_by_chr=args.output_bam_by_chr,
        tie_margin=args.tie_margin,
        no_realign=args.no_realign,
        no_polyA=args.no_polyA)

elif __name__ == '__main__':
    # Test units
    del sys.argv[1:] # Don't choke on extra command-line parameters
    import unittest
    import shutil
    import tempfile

    class TestPrintReadletizedOutput(unittest.TestCase):
        """ Tests print_readletized_output(). """
        def setUp(self):
            # For storing output of print_readletized_output()
            self.temp_dir_path = tempfile.mkdtemp()
            self.output_file = os.path.join(self.temp_dir_path,
                                                    'sample_output.tsv')

        def test_output(self):
            """ Fails if output of print_readletized_output() is wrong. """
            '''Every sequence has 100 bases.'''
            input_seqs = \
                    ['TTACATACCATACAGTGCGCTAGCGGGTGACAGATATAATGCAGATCCAT'
                     'ACAGACCAGATGGCAGACATGTGTTGCAGSCTGCAAGTGCAACGCGGTGA',
                     'GCAGAGTGCCGCAATGACGTGCGCCAAAGCGGTGACAGGGTGACAGTGAA'
                     'CCAAGTGACAAGTGAACAGGTGCCAGAGTGACCGAGTGACCAGTGGACCA',
                     'CAGAGTGCCGCAATGACGTGCGCCAAAGCGGACAAAGCACCATGACAAGT'
                     'ACACAGGTGACAGTGACAAGACAGAGGTGACACAGAGAAAGtGGGTGTGA',
                     'ATCGATTAAGCTATAACAGATAACATAGACATTGCGCCCATAATAGATAA'
                     'CTGACACCTGACCAGTGCCAGATGACCAGTGCCAGATGGACGACAGTAGC',
                     'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
                     'CTGACACCTGACCAGTGCCAGATGACCAGTGCCAGATGGACGACAGTAGC',
                     'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
                     'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA']
            sample_indexes = [
                    {'1' : 2}, {'1' : 2, '2' : 3}, {'1' : 4, '3' : 5},
                    {'1' : 5}, {'1' : 5}, {'1' : 5}
                ]
            reversed_complement_sample_indexes = [
                    {'1' : 7}, {'1' : 2, '2' : 2, '3' : 3},
                    {'1' : 8, '2' : 1, '3' : 6}, {'1' : 1},
                    {'1' : 1}, {'1' : 1}
                ]
            cap_sizes = []
            cap_size = 25
            while cap_size <= 50:
                cap_sizes.append(cap_size)
                cap_size = int(cap_size*2)
                if cap_size == cap_sizes[-1]:
                  cap_size += 1
            if cap_size != 50:
                # Always have a start or end read of length max_readlet_size
                cap_sizes.append(50)
            with open(self.output_file, 'w') as output_stream:
                for i in xrange(6):
                    print_readletized_output(
                            input_seqs[i], sample_indexes[i],
                            reversed_complement_sample_indexes[i],
                            '0:' + str(i), cap_sizes,
                            output_stream=output_stream,
                            min_readlet_size=25, readlet_interval=5,
                            max_readlet_size=50, no_polyA=True
                        )
            collected_readlets = []
            with open(self.output_file) as processed_stream:
                for readlet in processed_stream:
                    collected_readlets.append(readlet.rstrip().split('\t')[1:])
            '''Each read from input_reads spans 100 bases, and from the
            arguments passed to go() above, noncapping readlets should span 50.
            The capping fraction should arrange for one readlet spanning 25
            bases on either end of a given read. There should be 13 readlets in
            total. Spot-check some readlets after checking for read info.'''
            read_info = [info.split('\x1e') for _, info
                            in collected_readlets]
            read_info = [info[3:4]
                            + [info[-i].split('\x1f')
                                for i in xrange(4, 0, -1)]
                            for info in read_info if len(info) > 3]
            sample_indexes, reversed_complement_sample_indexes \
                = [{} for i in xrange(6)], [{} for i in xrange(6)]
            # All-A read should have been skipped completely
            self.assertEquals(len(read_info), 5)
            for i in xrange(5):
                for j in xrange(len(read_info[i][1])):
                    sample_indexes[i][read_info[i][1][j]] \
                        = int(read_info[i][3][j])
                    reversed_complement_sample_indexes[i][read_info[i][2][j]] \
                        = int(read_info[i][4][j])
                read_info[i] = read_info[i][:-2]
                read_info[i][1] = sample_indexes[i]
                read_info[i][2] = reversed_complement_sample_indexes[i]
            self.assertTrue([
                    'TTACATACCATACAGTGCGCTAGCGGGTGACAGATATAATGCAGATCCAT'
                    'ACAGACCAGATGGCAGACATGTGTTGCAGSCTGCAAGTGCAACGCGGTGA',
                    sample_indexes[0],
                    reversed_complement_sample_indexes[0]
                ] in read_info)
            self.assertTrue([
                    'GCAGAGTGCCGCAATGACGTGCGCCAAAGCGGTGACAGGGTGACAGTGAA'
                    'CCAAGTGACAAGTGAACAGGTGCCAGAGTGACCGAGTGACCAGTGGACCA',
                    sample_indexes[1],
                    reversed_complement_sample_indexes[1]
                ] in read_info)
            self.assertTrue([
                    'CAGAGTGCCGCAATGACGTGCGCCAAAGCGGACAAAGCACCATGACAAGT'
                    'ACACAGGTGACAGTGACAAGACAGAGGTGACACAGAGAAAGtGGGTGTGA',
                    sample_indexes[2],
                    reversed_complement_sample_indexes[2]
                ] in read_info)
            self.assertTrue([
                    'ATCGATTAAGCTATAACAGATAACATAGACATTGCGCCCATAATAGATAA'
                    'CTGACACCTGACCAGTGCCAGATGACCAGTGCCAGATGGACGACAGTAGC',
                    sample_indexes[3],
                    reversed_complement_sample_indexes[3]
                ] in read_info)
            self.assertTrue([
                    'ATCGATTAAGCTATAACAGATAACATAGACATTGCGCCCATAATAGATAA'
                    'CTGACACCTGACCAGTGCCAGATGACCAGTGCCAGATGGACGACAGTAGC',
                    sample_indexes[3],
                    reversed_complement_sample_indexes[3]
                ] in read_info)
            self.assertTrue([
                    'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
                    'CTGACACCTGACCAGTGCCAGATGACCAGTGCCAGATGGACGACAGTAGC',
                    sample_indexes[3],
                    reversed_complement_sample_indexes[3]
                ] in read_info)
            # Capping readlets
            self.assertTrue([       'ATCGATTAAGCTATAACAGATAACA'
                                    'TAGACATTGCGCCCATAATAGATAA',
                                    '0:3+\x1e0\x1e50'
                                ] in collected_readlets or
                                [   'ATCGATTAAGCTATAACAGATAACA'
                                    'TAGACATTGCGCCCATAATAGATAA',
                                    '0:3+\x1e0\x1e50\x1e'
                                    'ATCGATTAAGCTATAACAGATAACA'
                                    'TAGACATTGCGCCCATAATAGATAA'
                                    'CTGACACCTGACCAGTGCCAGATGA'
                                    'CCAGTGCCAGATGGACGACAGTAGC'
                                    '\x1e1\x1e1\x1e1\x1e1'
                                ] in collected_readlets
                            )
            self.assertTrue([       'CCAGTGCCAGATGGACGACAGTAGC',
                                    '0:3+\x1e75\x1e0'
                                ] in collected_readlets or
                                [   'CCAGTGCCAGATGGACGACAGTAGC',
                                    '0:3+\x1e75\x1e0\x1e'
                                    'ATCGATTAAGCTATAACAGATAACA'
                                    'TAGACATTGCGCCCATAATAGATAA'
                                    'CTGACACCTGACCAGTGCCAGATGA'
                                    'CCAGTGCCAGATGGACGACAGTAGC'
                                    '\x1e1\x1e1\x1e1\x1e1'
                                ] in collected_readlets
                            )
            # Noncapping readlets
            self.assertTrue([       'GTCAG'
                                    'TTATCTATTATGGGCGCAATGTCTA'
                                    'TGTTATCTGTTATAGCTTAA',
                                    '0:3-\x1e5\x1e45'
                                ] in collected_readlets or
                                [   'GTCAG'
                                    'TTATCTATTATGGGCGCAATGTCTA'
                                    'TGTTATCTGTTATAGCTTAA',
                                    '0:3-\x1e5\x1e45\x1e'
                                    'ATCGATTAAGCTATAACAGATAACA'
                                    'TAGACATTGCGCCCATAATAGATAA'
                                    'CTGACACCTGACCAGTGCCAGATGA'
                                    'CCAGTGCCAGATGGACGACAGTAGC'
                                    '\x1e1\x1e1\x1e1\x1e1'
                                ] in collected_readlets
                            )
            self.assertTrue([       'TAATAGATAA'
                                    'CTGACACCTGACCAGTGCCAGATGA'
                                    'CCAGTGCCAGATGGA',
                                    '0:3+\x1e40\x1e10'
                                ] in collected_readlets or
                                [   'TAATAGATAA'
                                    'CTGACACCTGACCAGTGCCAGATGA'
                                    'CCAGTGCCAGATGGA',
                                    '0:3+\x1e40\x1e10\x1e'
                                    'ATCGATTAAGCTATAACAGATAACA'
                                    'TAGACATTGCGCCCATAATAGATAA'
                                    'CTGACACCTGACCAGTGCCAGATGA'
                                    'CCAGTGCCAGATGGACGACAGTAGC'
                                    '\x1e1\x1e1\x1e1\x1e1'
                                ] in collected_readlets
                            )
            # Ensure no polyA readlets
            self.assertEquals([readlet for readlet in collected_readlets
                                if set(readlet[0]) == _polyA], [])

        def tearDown(self):
            # Kill temporary directory
            shutil.rmtree(self.temp_dir_path)

    unittest.main()
