#!/usr/bin/env python
"""
generate_bioreps.py
Abhi Nellore / July 14, 2014

Generates 20 bioreps for Rail paper using Flux Simulator by:
1) starting a single Flux simulation and interrupting the pipeline after
    6 columns of a PRO file P have been written.
2) reading RPKMs from Geuvadis data for by default a group a 20 YRI samples.
3) writing 20 files, each a copy of P corresponding to a different sample with
    a different coverage distribution. Column 6 has the absolute number of
    RNA molecules associated with each transcript. To compute this value,
    an RPKM for a given transcript from (1) is multiplied by 7. This
    gives ~10^6-10^7 nb_molecules for a given Flux simulation, which
    is close to that for the human example on the Flux website and apparently
    makes for adequate library yield.
4) restarts Flux, now for each PRO file, generating bioreps. FASTAs are named
    for samples.

The default locations of files from command-line parameters are on the Hopkins
Homewood High-Performance Cluster. Command-line parameter help specifies where
on the web to grab them and what was used for the paper.
"""
import argparse
from collections import defaultdict
import pandas as pd

if __name__ == '__main__':
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('-g', '--gtf', type=str,
        default=('/scratch0/langmead-fs1/geuvadis_sim/'
                 'gencode.v12.annotation.gtf')
        help=('Transcript annotation to pass to Flux; for Rail paper, this is '
              'the Gencode v12 GTF/GFF3 file obtainable at '
              'ftp://ftp.sanger.ac.uk/pub/gencode/release_12/'
              'gencode.v12.annotation.gtf.gz')
        )
    parser.add_argument('-f', '--flux', type=str,
        default=('/scratch0/langmead-fs1/shared/flux-simulator-1.2.1/bin/'
                 'flux-simulator')
        help=('Flux Simulator executable. v1.2.1 is used for paper and is '
              'obtainable at '
              'http://sammeth.net/artifactory/barna/barna/barna.simulator/'
              '1.2.1/flux-simulator-1.2.1.tgz')
        )
    parser.add_argument('-r', '--rpkm', type=str,
        default='/scratch0/langmead-fs1/geuvadis_sim/GD660.TrQuantRPKM.txt',
        help=('RPKM file with transcript labels (rows) and sample labels '
              '(columns). File is obtainable at '
              'http://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-3/files/'
              'analysis_results/GD660.TrQuantRPKM.txt')
        )
    parser.add_argument('-p', '--num-processes', type=int,
        default=5,
        help=('Number of instances of Flux Simulator to run simultaneously; '
              'set this to below 10 so it doesn\'t choke.')
        )
    parser.add_argument('-s', '--samples', type=str,
        default=('NA18508.1.M_111124_1,NA18510.3.M_120202_7,'
                 'NA18511.1.M_120209_1,NA18517.1.M_120209_7,'
                 'NA18519.1.M_111124_3,NA18520.1.M_111124_1,'
                 'NA18858.1.M_120209_7,NA18861.1.M_120209_2,'
                 'NA18907.1.M_120209_1,NA18908.7.M_120219_7,'
                 'NA18909.1.M_120209_8,NA18910.1.M_120209_5,'
                 'NA18912.1.M_120209_5,NA18916.1.M_120209_7,'
                 'NA18917.1.M_111124_5,NA18923.2.M_111216_5,'
                 'NA18933.1.M_111124_2,NA18934.1.M_120209_3,'
                 'NA19092.1.M_111124_3,NA19093.7.M_120219_7')
        help=('Comma-separated list of sample names from RPKM file whose '
              'expression profiles are to be mimicked. Defaults to list of '
              'sample names used in Rail paper (some YRIs) and determines '
              'number of simulations to perform. Be sure to exclude samples '
              'mentioned at http://geuvadiswiki.crg.es/index.php/'
              'QC_sample_info. Sample information is available at '
              'http://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/'
              'samples/.')
        )
    parser.add_argument('-c', '--read-count', type=int,
            default=40000000,
            help='Number of reads to generate per sample')
    parser.add_argument('--single-end', action='store_const', const=True,
            default=False,
            help='Generate single-end rather than paired-end reads')

    args = parser.parse_args()

    pd.DataFrame.from_csv(args.rpkm, sep=',')
    with open(args.rpkm) as rpkm_stream:
        sample_list = rpkm_stream.readline().strip().split('\t')[4:]
        samples = defaultdict(int)
        relevant_samples = set(args.samples.strip().split(','))
        # Store column index of each relevant sample
        for i, sample in enumerate(sample_list):
            if sample in relevant_samples:
                samples[sample] = i + 4
        # For storing row indexes of transcripts
        transcript_to_row = defaultdict(int)
        for i, line in enumerate(rpkm_stream):
