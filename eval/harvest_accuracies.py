#!/usr/bin/env python
"""
harvest_accuracies.py

Constructs tables of accuracy results from the files
performance_summary, performance_mapping_accuracy_summary,
performance_mapping_accuracy_SC_summary, and
performance_intron_recovery_summary contained in output directories of
of either grab_112_sim_results.sh or the commands output by
create_commands_for_all_sample_sims.sh.

A tab-separated matrix A_ij is output for each of intron recovery performance,
exon-exon junction recovery performance, mapping accuracy performance,
and mapping accuracy performance where all alignments for which > 10 percent of
bases are soft-clipped are treated as unmapped. Each i is a different mode
(of a given aligner), and each j is a sample name, mean, or stdev.
A_ij is the comma-separated list (precision, recall, f-score)
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

def stats(path, first=True):
    """ Returns comma-separated list of accuracy stats
        for summary files output by intron_recovery_performance.py,
        spliced_read_recovery_performance.py, and mapping_accuracy.py.

        Included are (precision, recall, f-score)
        
        path: path to perform file
        first: whether to choose first or second precision/recall from
            performance file. Relevant when reading output of
            mapping_accuracy.py, for which the first value applies to
            basewise performance metrics while the second applies to
            read-level performance metrics

        
        Return value: (precision, recall, f-score)
    """
    if first:
        index = 1
    else:
        index = 2
    with open(path) as accuracy_stream:
        for line in accuracy_stream:
            tokens = line.strip().split('\t')
            if tokens[0] == 'precision':
                precision = float(tokens[index])
            elif tokens[0] == 'recall':
                recall = float(tokens[index])
    fscore = 2 * precision * recall / (precision + recall)
    return (precision, recall, fscore)

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

def write_samples(samples):
    """ Writes list of samples to stdout for matrix

        samples: list of sample names

        No return value
    """
    print '\t'.join([''] + samples + ['means', 'stdevs'])

if __name__ == '__main__':
    import argparse
    import subprocess
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
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
    samples = args.samples.split('\t')
    for perform in ['perform_intron_recovery_summary', 'perform_summary']:
        print perform + '\n'
        write_samples(samples)
        for mode in modes:
            sample_stats = []
            for sample in samples:
                full_path = os.path.join(args.sam_dir, sample, mode, perform)
                sample_stats.append(stats(full_path))
            sample_stats = zip(*sample_stats)
            means = [mean(stat_list) for stat_list in sample_stats]
            stdevs = [stdev(stat_list) for stat_list in sample_stats]
            sample_stats = zip(*sample_stats)
            sys.stdout.write(mode)
            for i, sample in enumerate(samples):
                sys.stdout.write('\t' + ('%.3f,%.3f,%.3f' % sample_stats[i]))
            sys.stdout.write('\n')
    for perform in ['perform_mapping_accuracy_summary',
                    'perform_mapping_accuracy_SC_summary']:
        print perform + '\n'
        for first in [False, True]:
            if first:
                print 'basewise'
            else:
                print 'readwise'
            write_samples(samples)
            for mode in modes:
                sample_stats = []
                for sample in samples:
                    full_path = os.path.join(
                                        args.sam_dir, sample, mode, perform
                                    )
                    sample_stats.append(stats(full_path), first=first)
                sample_stats = zip(*sample_stats)
                means = [mean(stat_list) for stat_list in sample_stats]
                stdevs = [stdev(stat_list) for stat_list in sample_stats]
                sample_stats = zip(*sample_stats)
                sys.stdout.write(mode)
                for i, sample in enumerate(samples):
                    sys.stdout.write(
                                '\t' + ('%.3f,%.3f,%.3f' % sample_stats[i])
                            )
                sys.stdout.write('\n')
