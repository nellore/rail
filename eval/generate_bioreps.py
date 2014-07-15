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
    and subsequently multiplied by 10. This gives ~1e6-1e7 nb_molecules for a
    given Flux simulation, which is close to that for the human example on the
    Flux website and apparently makes for adequate library yield.
4) restarting Flux, now for each PRO file, generating bioreps. FASTAs are named
    for samples.

The default locations of files from command-line parameters are on the Hopkins
Homewood High-Performance Cluster. Command-line parameter help specifies where
on the web to grab them and what parameters are used in the paper.

REQUIRES PANDAS. Download the Anaconda distribution of Python to simplify
getting it.

NOTE that expressed fractions (column 5 of pro files) are incorrect, but they
are never used by Flux Simulator. Also note that the script "assumes" a direct
proportionality between RPKMs*read length and number of RNA molecules -- this
isn't even true in simulation. In fact, the simulated correlation appears to be
between 80 and 90 percent. (To reproduce this, compute the correlation between
the 6th and 10th columns of a given pro file. awk code:

cat X.pro | awk '{print $6 "\t" $10}' 
| awk '{ xy+=($1*$2); x+=$1; y+=$2; x2+=($1*$1); y2+=($2*$2); } 
END { print "NR=" NR; ssx=x2-((x*x)/NR); 
print "ssx=" ssx; ssy=y2-((y*y)/NR); print "ssy=" ssy; 
ssxy = xy - ((x*y)/NR); print "ssxy=" ssxy; r=ssxy/sqrt(ssx*ssy); 
print "r=" r; }'.)

based on http://awk.info/?doc/tools/correlate.html .)
The intent here, however, is to create a relatively
realistic coverage distribution and capture something like the variation
observed in what users might regard as a set of bioreplicates.

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

def write_par_and_pro(par_template, pro_template, basename,
                      seed, rpkms, sample):
    """ Writes PAR and PRO file for a given sample.

        args: ordered tuple whose components are:
        par_template: PAR template file
        pro_template: PRO template file
        basename: Destination PAR/PRO basename
        seed: Random seed to put in PAR (integer)
        rpkms: Pandas dataframe encoding transcript RPKMs across samples
        sample: Sample name

        Return value: 0.
    """
    with open(basename + '.pro', 'w') as write_stream:
        with open(pro_template) as read_stream:
            for line in read_stream:
                tokens = line.strip().split('\t')
                transcript_label = tokens[1]
                transcript_length_in_kb = float(tokens[3]) / 1000
                try:
                    rpkm = rpkms.loc[transcript_label][sample]
                except KeyError:
                    # Not expressed; kill it
                    rpkm = 0.0
                tokens[5] = transcript_length_in_kb * rpkm * 10
                remaining = tokens[5] - int(tokens[5])
                if remaining > 0.5:
                    tokens[5] = int(tokens[5]) + 1
                else:
                    tokens[5] = int(tokens[5])
                tokens[5] = str(tokens[5])
                print >>write_stream, '\t'.join(tokens)
    # Write distinct seed for sample
    with open(basename + '.par', 'w') as write_stream:
        with open(par_template) as read_stream:
            # Next-to-last param is REF_FILE_NAME; last is SEED
            all_parameters = read_stream.read().strip().split('\n')
            ref_file_name = all_parameters[-2].split('\t')
            assert ref_file_name[0] == 'REF_FILE_NAME'
            sorted_name = ref_file_name[1][:-4] + '_sorted.gtf'
            if os.path.exists(sorted_name):
                name_to_write = sorted_name
            else:
                name_to_write = ref_file_name[1]
            print >>write_stream, \
                '\n'.join(all_parameters[:-2])
            print >>write_stream, 'REF_FILE_NAME\t%s' % name_to_write
            print >>write_stream, 'SEED\t%d' % seed
    return 0

def run_flux(par, flux):
    """ Runs Flux pipeline after creation of PRO file.

        par: PAR file; there must be a corresponding PRO file in the
            same directory
        flux: flux executable

        Return value: Flux exitlevel.
    """
    with open(par + '.log', 'w') as log_stream:
        return_value = subprocess.call([flux, '-l', '-s', '-p', par],
                                        stderr=log_stream,
                                        stdout=open(os.devnull, 'w'))
    return return_value

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
    	default='./',
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
        default=7,
        help=('Number of instances of Flux Simulator to run simultaneously; '
              'set this to below 10 so it doesn\'t choke.')
        )
    parser.add_argument('--fasta', type=str,
        default='/scratch0/langmead-fs1/shared/references/hg19/fasta',
        help='Where to find reference FASTAs for chrs. Use GRCh37.')
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

    print >>sys.stderr, 'Reading RPKMs...'
    rpkms = pd.DataFrame.from_csv(args.rpkm, sep='\t')
    relevant_samples = args.samples.strip().split(',')

    temp_dir = tempfile.mkdtemp()
    atexit.register(kill_dir, temp_dir)
    expression_par = os.path.join(temp_dir, 'sim.par')
    par_template = [
        ('NB_MOLECULES', '5000000'),
        ('LOAD_NONCODING', 'YES'),
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
        ('UNIQUE_IDS', 'YES'),
        ('GEN_DIR', args.fasta),
        ('REF_FILE_NAME', args.gtf)]
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
    for i, sample in enumerate(relevant_samples):
        pool.apply_async(write_par_and_pro,
                       (expression_par, expression_pro,
                            os.path.join(args.output, sample + '_sim'),
                            i, rpkms, sample),
                       callback=return_values.append)
    pool.close()
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
    for sample in relevant_samples:
        pool.apply_async(run_flux,
                         (os.path.join(args.output, sample + '_sim.par'),
                          args.flux),
                         callback=return_values.append)
    pool.close()
    while len(return_values) != relevant_count:
        sys.stdout.write('Completed %d/%d sims.\r'
                            % (len(return_values), relevant_count))
        sys.stdout.flush()
        time.sleep(.2)
    print >>sys.stderr, 'Done.'