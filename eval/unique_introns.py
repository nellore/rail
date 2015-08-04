#!/usr/bin/env python
"""
unique_introns.py

Measures precision and recall of introns that are unique to each of the 20
samples selected at random (in create_single_sample_sim_commands.py) from
112 simulated GEUVADIS-based samples. Assumes a directory structure for output
created by either grab_112_sim_results.sh or the commands output by
create_commands_for_all_sample_sims.sh.

A file Aligned.out.sam is assumed to be in output directories for all aligners
but Rail-RNA and TopHat 2. TopHat 2's output file is assumed to be
accepted_hits.bam, and Rail's output is assumed to span all the bam files in
its output directory.

A tab-separated matrix A_ij is output. Each i is a different mode
(of a given aligner), and each j is a sample name, mean, or stdev.
A_ij is the comma-separated list (true intron count, retrieved intron count,
                                    overlap, precision, recall, f-score)
"""

import os
from collections import defaultdict
import sys
import glob
import time
import math
import multiprocessing
import signal
import re
from intron_recovery_performance import indels_introns_and_exons

def init_worker():
    """ Prevents KeyboardInterrupt from reaching a pool's workers.

        Exiting gracefully after KeyboardInterrupt or SystemExit is a
        challenge. The solution implemented here is by John Reese and is from
        http://noswap.com/blog/python-multiprocessing-keyboardinterrupt .

        No return value.
    """
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def introns_from_bed(index_and_bed):
    """ Converts BED to dictionary that maps RNAMES to sets of introns.
        
        index_and_bed is composed of:
        index: index of bed file; used for tracking which bed is which
            when using multiple cores
        bed: input BED file characterizing splice junctions

        Return value: (index, a dictionary). Each key is an RNAME, typically a
            chromosome, and its corresponding value is a set of tuples, each
            denoting an intron on RNAME. Each tuple is of the form
            (start position, end position).
    """
    index, bed = index_and_bed
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
    return (index, introns)

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

def introns_from_sam(index_and_sam_and_samtools_exe):
    """ Writes output that maps QNAMES to exon-exon junctions overlapped.

        index_and_sam_and_samtools_exe is composed of
        index: index of SAM file; used for tracking which bed is which
            when using multiple cores
        sam: where to find retrieved alignments in SAM form
        samtools_exe: where to find samtools executable

        Return value: (index, a dictionary). Each key is an RNAME, typically a
            chromosome, and its corresponding value is a set of tuples, each
            denoting an intron on RNAME. Each tuple is of the form
            (start position, end position).
    """
    index, sam, samtools_exe = index_and_sam_and_samtools_exe
    introns = set()
    for sam_file in sam:
        sam_stream = subprocess.Popen([samtools_exe, 'view', sam_file],
                                            stdout=subprocess.PIPE)
        for line in sam_stream.stdout:
            if line[0] == '@': continue
            try:
                tokens = line.strip().split('\t')
                flag = int(tokens[1])
                if flag & 4:
                    continue
                name = tokens[0]
                rname = tokens[2]
                cigar = tokens[5]
                pos = int(tokens[3])
                seq = tokens[9]
                flag = int(tokens[1])
                if 'N' not in cigar or flag & 256:
                    continue
                _, _, introns_to_add, _ = indels_introns_and_exons(cigar,
                                            dummy_md_index(cigar), pos, seq)
                for intron in introns_to_add:
                    introns.add((rname,) + intron[:2])
            except IndexError:
                print >>sys.stderr, ('Error found on line: ' + line)
                raise
            except Exception as e:
                print >>sys.stderr, e
                raise
        sam_stream.kill()
    return (index, introns)

def stats(index_and_true_and_retrieved):
    """ Returns comma-separated list of accuracy stats

        Included are (relevant instances, retrieved instances, overlap,
                        precision, recall, f-score)
        
        index_and_true_and_retrieved is composed of:
        index: used to track which sample is associated with the stats
        true: set of true instances
        retrieved: set of retrieved indexes

        Return value: (index, tuple of stats)
    """
    index, true, recovered = index_and_true_and_retrieved
    relevant = len(true)
    retrieved = len(recovered)
    overlap = len(true.intersection(recovered))
    try:
        precision = float(overlap) / retrieved
    except ZeroDivisionError:
        precision = 1
    try:
        recall = float(overlap) / relevant
    except ZeroDivisionError:
        recall = 1
    try:
        fscore = 2 * precision * recall / (precision + recall)
    except ZeroDivisionError:
        fscore = 1
    return (index, (relevant, retrieved, overlap, precision, recall, fscore))

def mean(a_list):
    """ Finds mean of list.

        a_list: a list

        Return value: mean
    """
    return float(sum(a_list)) / len(a_list)

def stdev(a_list):
    """ Finds sample stdev of list.

        a_list: a list

        Return value: stdev
    """
    working_mean = mean(a_list)
    return math.sqrt(1. / (len(a_list) - 1)
                        * sum((el-working_mean)**2 for el in a_list))

if __name__ == '__main__':
    import argparse
    import subprocess
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--samtools', type=str, required=False,
            default='samtools',
            help=('Path to SAMTools executable; used to extract introns from '
                  'BAM files')
        )
    parser.add_argument('--samples', type=str, required=False,
            default=('HG00115_male_GBR_UU_6-1-1,'
                     'NA06984_male_CEU_UNIGE_1-1-1,'
                     'HG00313_female_FIN_UNIGE_1-1-1,'
                     'NA19095_female_YRI_LUMC_7-1-3,'
                     'HG00117_male_GBR_LUMC_7-1-1,'
                     'NA20768_female_TSI_HMGU_5-1-1,'
                     'NA20582_female_TSI_ICMB_4-1-1,'
                     'NA19130_male_YRI_HMGU_5-1-1,'
                     'HG00139_male_GBR_LUMC_7-1-1,'
                     'NA18486_male_YRI_LUMC_7-1-1,'
                     'HG01790_female_GBR_MPIMG_3-1-1,'
                     'NA12287_female_CEU_UNIGE_1-1-1,'
                     'NA12287_female_CEU_MPIMG_3-1-2,'
                     'HG00096_male_GBR_UNIGE_1-1-1,'
                     'NA11831_male_CEU_LUMC_7-1-1,'
                     'NA12874_male_CEU_UNIGE_1-1-1,'
                     'HG00154_female_GBR_HMGU_5-1-1,'
                     'NA07051_male_CEU_HMGU_5-1-1,'
                     'NA12776_female_CEU_UU_6-1-1,'
                     'NA20813_female_TSI_HMGU_5-1-1'),
            help=('Comma-separated list of sample names to use in '
                  'computation of introns unique to particular samples.')
        )
    parser.add_argument('-t', '--true-introns-bed-dir', type=str,
            required=True,
            help='Path to directory containing Flux Simulator BEDs'
        )
    parser.add_argument('-r', '--sam-dir', type=str,
            required=False,
            default=None,
            help=('Path to root directory containing sample dirs with aligner '
                  'results')
        )
    parser.add_argument('-a', '--aligners', type=str,
            required=True, nargs='+',
            help=('Aligners whose results this script should look for. Use '
                  'hisat, rail, subjunc, tophat, or star.')
        )
    args = parser.parse_args()
    # Ingest true introns
    true_introns = defaultdict(list)
    samples = args.samples.split(',')
    bed_paths = [os.path.join(args.true_introns_bed_dir,
                                        sample + '_sim.bed')
                            for sample in samples]
    print >>sys.stderr, 'Loading true introns...'
    # Use multiple cores
    pool = multiprocessing.Pool(multiprocessing.cpu_count() - 1,
                                    init_worker, maxtasksperchild=5)
    returned_introns = []
    pool.map_async(
                    introns_from_bed,
                    list(enumerate(bed_paths)),
                    callback=returned_introns.extend
                )
    bed_path_size = len(bed_paths)
    while len(returned_introns) < bed_path_size:
        time.sleep(1)
    for index, true_intron_set in returned_introns:
        for intron in true_intron_set:
            true_introns[intron].append(index)
    true_unique_introns = defaultdict(set)
    for intron in true_introns:
        if len(true_introns[intron]) == 1:
            true_unique_introns[true_introns[intron][0]].add(intron)
    print >>sys.stderr, 'Loaded true introns.'
    print '\t'.join([''] + samples + ['means', 'stdevs'])
    modes = []
    if 'star' in args.aligners:
        modes.extend(['star/ann_paired_1pass',
                        'star/ann_paired_2pass',
                        'star/noann_paired_1pass',
                        'star/nogen_noann_paired_2pass'])
    if 'tophat' in args.aligners:
        modes.extend(['tophat/ann_paired', 'tophat/noann_paired'])
    if 'subjunc' in args.aligners:
        modes.extend(['subjunc'])
    if 'rail' in args.aligners:
        modes.extend(['rail'])
    if 'hisat' in args.aligners:
        modes.extend(['hisat/ann_paired_1pass',
                        'hisat/ann_paired_2pass',
                        'hisat/noann_paired_1pass',
                        'hisat/noann_paired_2pass'])
    for mode in modes:
        sys.stdout.write(mode)
        if 'rail' in mode:
            alignment_files = [glob.glob(
                    os.path.join(args.sam_dir,
                                    sample, mode,
                                    'alignments/alignments.*.bam')
                ) for sample in samples]
            if not alignment_files[0]:
                # working with 112 sim
                alignment_files = [glob.glob(
                            os.path.join(args.sam_dir,
                                            sample, 'alignments.*.bam')
                        ) for sample in samples]
        elif 'tophat' in mode:
            alignment_files = [
                    [os.path.join(args.sam_dir, sample,
                                    mode, 'accepted_hits.bam')]
                    for sample in samples
                ]
        else:
            alignment_files = [
                    [os.path.join(args.sam_dir, sample,
                                    mode, 'Aligned.out.sam')]
                    for sample in samples
                ]
        returned_introns = []
        retrieved_introns = defaultdict(list)
        pool.map_async(
                    introns_from_sam,
                    [el + (args.samtools,) for el
                        in enumerate(alignment_files)],
                    callback=returned_introns.extend
                )
        alignment_files_size = len(alignment_files)
        while len(returned_introns) < alignment_files_size:
            time.sleep(1)
        for index, retrieved_intron_set in returned_introns:
            for intron in retrieved_intron_set:
                retrieved_introns[intron].append(index)
        retrieved_unique_introns = defaultdict(set)
        for intron in retrieved_introns:
            if len(retrieved_introns[intron]) == 1:
                retrieved_unique_introns[
                        retrieved_introns[intron][0]
                    ].add(intron)
        returned_stats = []
        pool.map_async(
                        stats,
                        [(i, true_unique_introns[i],
                          retrieved_unique_introns[i])
                         for i in xrange(len(samples))],
                        callback=returned_stats.extend
                    )
        sample_size = len(samples)
        while len(returned_stats) < sample_size:
            time.sleep(1)
        means = tuple([mean([el[1][k] for el in returned_stats])
                       for k in xrange(6)])
        stdevs = tuple([
                stdev([el[1][k] for el in returned_stats]) for k in xrange(6)
            ])
        returned_stats.sort(key=lambda x: x[0])
        for stat_set in returned_stats:
            sys.stdout.write('\t%d,%d,%d,%.3f,%.3f,%.3f' % stat_set[1])
        sys.stdout.write('\t%.3f,%.3f,%.3f,%.3f,%.3f,%.3f' % means)
        sys.stdout.write('\t%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n' % stdevs)
        sys.stdout.flush()
    pool.close()
    pool.join()
