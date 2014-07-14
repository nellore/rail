#!/usr/bin/env python
"""
generate_bioreps.py
Abhi Nellore / July 14, 2014

Generates 20 bioreps for Rail paper using Flux Simulator by:

1) starting a single Flux simulation and interrupting the pipeline after
    6 columns of a PRO file P have been written.
2) reading RPKMs from Geuvadis data for by default a group a 20 YRI samples.
3) writing 20 files, each a copy of P except with a different coverage
    distribution. The coverage distribution is specified as follows: column 6
    of P has the absolute number of RNA molecules associated with each
    transcript. To recompute this value, an RPKM for a given transcript from
    (2) is multiplied by transcript length (column 4 of P) in kb
    and subsequently doubled. This gives ~10^6-10^7 nb_molecules for a given
    Flux simulation, which is close to that for the human example on the Flux
    website and apparently makes for adequate library yield.
4) restarting Flux, now for each PRO file, generating bioreps. FASTAs are named
    for samples.

The default locations of files from command-line parameters are on the Hopkins
Homewood High-Performance Cluster. Command-line parameter help specifies where
on the web to grab them and what parameters are used in the paper.

REQUIRES PANDAS. Download the Anaconda distribution of Python to simplify
getting it.

Run this file from desired output directory if not entering output dir with -o.
"""

import os
import sys
import argparse
from collections import defaultdict
import pandas as pd
import subprocess
import tempfile
import multiprocessing
import time
import atexit

def kill_dir(path_to_dir):
    """ Removes directory tree.

        path_to_dir: directory to remove

        No return value.
    """
    try:
        shutil.rmtree(path_to_dir)
    except: pass

def write_par_and_pro(args):
    """ Writes PAR and PRO file for a given sample.

        args: ordered tuple whose components are:
                1) PAR template file
                2) PRO template file
                3) Destination PAR/PRO basename
                4) Random seed to put in PAR (integer)
                5) Pandas dataframe encoding transcript RPKMs across samples
                6) Sample name

        Return value: 0.
    """
    par_template, pro_template, basename, seed, rpkms, sample = args
    with open(basename + '.pro', 'w') as write_stream:
        with open(pro_template) as read_stream:
            for line in read_stream:
                tokens = line.strip().split('\t')
                transcript_label = tokens[0]
                transcript_length_in_kb = float(tokens[4]) / 1000
                try:
                    rpkm = rpkms.loc(transcript_label)[sample]
                except KeyError:
                    # Not expressed; kill it
                    rpkm = 0
                tokens[5] = str(transcript_length_in_kb * rpkm * 2)
                print >>write_stream, '\t'.join(tokens)
    # Write distinct seed for sample
    with open(basename + '.par', 'w') as write_stream:
        with open(par_template) as read_stream:
            print >>write_stream, read_stream.read()
            print >>write_stream, 'SEED\t%d' % seed
    return 0

def run_flux(args):
    """ Runs Flux pipeline after creation of PRO file.

        args: ordered tuple whose components are:
                par: PAR file; there must be a corresponding PRO file in the
                    same directory
                flux: flux executable

        Return value: Flux exitlevel.
    """
    par, flux = args
    with open(par + '.log', 'w') as log_stream:
        return subprocess.call([flux, '-l', '-s', '-p', par],
                                    stderr=log_stream,
                                    stdout=open(os.devnull, 'w'))

if __name__ == '__main__':
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('-g', '--gtf', type=str,
        default=('/scratch0/langmead-fs1/geuvadis_sim/'
                 'gencode.v12.annotation.gtf'),
        help=('Transcript annotation to pass to Flux; for Rail paper, this is '
              'the Gencode v12 GTF/GFF3 file obtainable at '
              'ftp://ftp.sanger.ac.uk/pub/gencode/release_12/'
              'gencode.v12.annotation.gtf.gz')
        )
    parser.add_argument('-o', '--output', type=str,
    	default='./'
        help='Where to put simulation FASTAs and BEDs')
    parser.add_argument('-f', '--flux', type=str,
        default=('/scratch0/langmead-fs1/shared/flux-simulator-1.2.1/bin/'
                 'flux-simulator'),
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
                 'NA19092.1.M_111124_3,NA19093.7.M_120219_7'),
        help=('Comma-separated list of sample names from RPKM file whose '
              'expression profiles are to be mimicked. Defaults to list of '
              'sample names used in Rail paper (some YRIs) and determines '
              'number of simulations to perform. Be sure to exclude samples '
              'mentioned at http://geuvadiswiki.crg.es/index.php/'
              'QC_sample_info. Sample information is available at '
              'http://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/'
              'samples/')
        )
    parser.add_argument('-c', '--read-count', type=int,
            default=40000000,
            help='Number of reads to generate per sample')
    parser.add_argument('--single-end', action='store_const', const=True,
            default=False,
            help='Generate single-end rather than paired-end reads')
    parser.add_argument('-l', '--read-length', type=int,
            default=76,
            help='Length of a generated read; WARNING: error model used by '
            	 'Flux here applies to reads 76 bases long')

    args = parser.parse_args()

    rpkms = pd.DataFrame.from_csv(args.rpkm, sep='\t')
    relevant_samples = args.samples.strip().split(',')

    temp_dir = tempfile.mkdtemp()
    atexit.register(kill_dir, temp_dir)
    expression_par = os.path.join(temp_dir, 'sim.par')
    par_template = [
        ('NB_MOLECULES', '5000000'),
        ('REF_FILE_NAME', args.gtf),
        ('GEN_DIR', temp_dir),
        ('LOAD_NONCODING', 'NO'),
        ('TSS_MEAN', '50'),
        ('POLYA_SCALE', 'NaN'),
        ('POLYA_SHAPE', 'NaN'),
        ('FRAG_SUBSTRATE', 'RNA'),
        ('FRAG_METHOD', 'UR'),
        ('FRAG_UR_ETA', '350'),
        ('FRAG_UR_D0', '1'),
        ('RTRANSCRIPTION', 'YES'),
        ('RT_PRIMER', 'RH'),
        ('RT_LOSSLESS', 'YES'),
        ('RT_MIN', '500'),
        ('RT_MAX', '5500'),
        ('GC_MEAN', 'NaN'),
        ('PCR_PROBABILITY', '0.050000'),
        ('FILTERING', 'NO'),
        ('READ_NUMBER', str(args.read_count)),
        ('READ_LENGTH', str(args.read_length)),
        ('PAIRED_END', 'NO' if args.single_end else 'YES'),
        ('ERR_FILE', '76'),
        ('FASTA', 'YES'),
        ('UNIQUE_IDS', 'YES')]
    with open(expression_par, 'w') as par_stream:
        print >>par_stream, '\n'.join(['\t'.join(parameter)
                                       for parameter in par_template])
    # Make reproducible by specifying random seed
    print >>par_stream, 'SEED\t0'
    print >>sys.stderr, 'Creating PRO template with Flux Simulator...'
    flux_expression_return = subprocess.call(
                                    [args.flux, '--force', '-x', '-p',
                                        expression_par],
                                    stderr=sys.stderr,
                                    stdout=sys.stdout
                                )
    if flux_expression_return:
        raise RuntimeError('Flux Simulator returned exitlevel %s.'
                            % flux_expression_return)
    expression_pro = expression_par[:-4] + '.pro'
    if not os.path.exists(expression_pro):
        raise RuntimeError('PRO template with same basename as PAR template '
                           'was not created.')
    print >>sys.stderr, 'Creating PAR and PRO files for bioreplicate sims...'
    pool = multiprocessing.Pool(args.num_processes)
    return_values = []
    pool.map_async(write_par_and_pro,
                    [(expression_par, expression_pro,
                        os.path.join(args.output, sample + '_sim'),
                        i, rpkms, sample) for i, sample
                        in enumerate(relevant_samples)],
                    callback=return_values.extend)
    relevant_count = len(relevant_samples)
    while len(return_values) != relevant_count:
        sys.stdout.write('Created %d/%d PAR/PRO pairs.\r' \
                            % (len(return_values), relevant_count))
        sys.stdout.flush()
        time.sleep(.2)
    print >>sys.stderr, 'Created all PAR/PRO pairs.'
    print >>sys.stderr, 'Running sims...'
    pool = multiprocessing.Pool(args.num_processes)
    return_values = []
    pool.map_async(run_flux,
                    [(sample + '_sim.par', args.flux)
                        for sample in relevant_samples])
    while len(return_values) != relevant_count:
        print 'Completed %d/%d sims.' \
            % (len(return_values), relevant_count)
        time.sleep(.2)
    print >>sys.stderr, 'Done.'