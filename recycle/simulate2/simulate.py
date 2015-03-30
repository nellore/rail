## quick RNA-seq DE simulator
## AF december 12, 2013
## BTL changes January 10, 2014

# TODO:
# - Generate multiple bio reps
# - Make simulated read names meaningful, pointing to point of origin
#   so we can analyze later
# - Kill rest of the hard-coded variables, e.g. baseline_expressions

## dependencies:
# pip install biopython
# pip install numpy

__author__ = "Alyssa Frazee and Ben Langmead"

from random import randint, uniform, choice, sample, gauss
try:
    import numpypy as np
except ImportError:
    pass
import numpy as np
from reference import ReferenceIndexed
from annotation import GeneAnnotation
import string
import logging

# TODO: these should be turned into command-line arguments as well

##### baseline expression levels for transcripts
##### e.g., what % of transcripts should be lowly expressed, medium-ly expressed, highly expressed?
fraction_high = 0.1
fraction_medium = 0.3
fraction_low = 0.1
almost_no_expression = 0.1
baseline_expressions = [100, 300, 900]  # # reads coming from low, medium, and high expression transcripts respectively

##### what fold changes should we see in the differentially expressed transcripts?
fold_change_opts = [2, 5]

# function for reverse-complementing
_revcomp_trans = string.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(s):
    return s[::-1].translate(_revcomp_trans)


def add_errors(s, nerrors):
    if nerrors == 0:
        return s
    ln = len(s)
    error_positions = sorted(sample(xrange(ln), nerrors))
    seq_list = []
    last_err = -1
    for error_position in error_positions:
        seq_list.append(s[last_err+1:error_position])
        new_nuc = s[error_position]
        while new_nuc == s[error_position]:
            new_nuc = choice('ACGT')
        seq_list.append(new_nuc)
        last_err = error_position
    seq_list.append(s[error_positions[-1]+1:])
    return ''.join(seq_list)


def simulate_from_transcript(tx, fragment_mean, fragment_sd, readlen, paired, error_rate, stranded=False):
    fraglen = int(round(gauss(fragment_mean, fragment_sd)))
    if tx.nucleotide_length() > fraglen:
        st = randint(0, tx.nucleotide_length() - fraglen)
        frag = tx.unspliced_substring_from_5p(st, fraglen)
        assert len(frag) == fraglen
    else:
        frag = tx.unspliced_substring_from_5p(0, tx.nucleotide_length())
        assert len(frag) <= fraglen
    m1_nerrors = np.random.binomial(readlen, error_rate)
    m1_error_seq = add_errors(frag[:readlen], m1_nerrors)
    m1_qual = 'I' * readlen
    if paired:
        m2_nerrors = np.random.binomial(readlen, error_rate)
        m2_error_seq = add_errors(reverse_complement(frag[-readlen:]), m2_nerrors)
        m2_qual = 'I' * readlen
    if not stranded and randint(0, 1) == 0:
        # Make fragment from anti-sense strand
        if paired:
            m1_error_seq, m2_error_seq = reverse_complement(m2_error_seq), reverse_complement(m1_error_seq)
        else:
            m1_error_seq = reverse_complement(m1_error_seq)
    name = 'simread'
    seq = m1_error_seq if not paired else (m1_error_seq, m2_error_seq)
    quals = m1_qual if not paired else (m1_qual, m2_qual)
    return name, seq, quals


def run_sequencer(transcripts, num_replicates, outfile_prefix, fragment_mean, fragment_sd, readlen, paired, error_rate, dispersion_param,
                  stranded=False):
    if paired:
        f1s = [open("%s_rep%d_1.fq" % (outfile_prefix, repi+1), 'w') for repi in xrange(num_replicates)]
        f2s = [open("%s_rep%d_2.fq" % (outfile_prefix, repi+1), 'w') for repi in xrange(num_replicates)]
    else:
        f1s = [open("%s_rep%d.fq" % (outfile_prefix, repi+1), 'w') for repi in xrange(num_replicates)]
        f2s = []
    nread_tot = 0
    for txi, tx_tup in enumerate(transcripts):
        tx, basemean, fold_change = tx_tup
        if (txi % 100) == 0:
            logging.info('  Generated %d reads from %d transcripts so far' % (nread_tot, txi))
        for repi in xrange(num_replicates):
            numreads = np.random.negative_binomial(n=dispersion_param,
                                                   p=float(dispersion_param) / (dispersion_param + basemean))
            nread_tot += numreads
            for read in xrange(numreads):
                name, seq, qual = simulate_from_transcript(tx, fragment_mean, fragment_sd, readlen, paired, error_rate,
                                                           stranded=stranded)
                if paired:
                    m1, m2 = seq
                    q1, q2 = qual
                    f1s[repi].write('@' + name + '/1' + '\n' + m1 + '\n+\n' + q1 + '\n')
                    f2s[repi].write('@' + name + '/2' + '\n' + m2 + '\n+\n' + q2 + '\n')
                else:
                    f1s[repi].write('@' + name + '\n' + seq + '\n+\n' + qual + '\n')
    for fh in f1s + f2s:
        fh.close()


def go(args):
    ###############################################
    ############## begin simulation ###############
    ###############################################

    logging.info('Parsing indexed FASTA...')
    ref = ReferenceIndexed(args.reference_fasta)

    logging.info('Parsing gene annotations...')
    annot = GeneAnnotation.from_gtf_file(args.gtf)

    logging.info('Populating gene annotations with sequences from genome...')
    annot.set_sequence(ref)

    logging.info('Summarizing gene annotations...')
    annot.summarize()

    info_file = open('simulation_information.txt', 'w')
    info_file.write('txid\tbase_expression\tDE\tfold_change\n')
    logging.info('Parsing transcripts...')
    txi = 0
    transcript_set = []
    for tx_id, tx_tup in annot.xscripts.iteritems():
        gene_id, tx = tx_tup
        meanflip = uniform(0, 1)
        if meanflip < fraction_low:
            basemean = baseline_expressions[0]
        elif meanflip < fraction_low + fraction_medium:
            basemean = baseline_expressions[1]
        elif meanflip < fraction_low + fraction_medium + fraction_high:
            basemean = baseline_expressions[2]
        else:
            basemean = almost_no_expression
        de = uniform(0, 1) < args.de_fraction
        fold_change = choice(fold_change_opts) if de else 1
        transcript_set.append((tx, basemean, fold_change))
        info_file.write(tx_id+'\t'+str(basemean)+'\t'+str(de)+'\t'+str(fold_change)+'\n')
        txi += 1
        if (txi % 5000) == 0:
            logging.info('  Assigned abundances to %d transcripts so far' % txi)
    info_file.close()

    # simulate reads!
    run_sequencer(transcripts=transcript_set, num_replicates=args.bio_replicates, outfile_prefix=args.output_prefix,
                  fragment_mean=args.fragment_mean, fragment_sd=args.fragment_sd, readlen=args.read_length,
                  paired=args.paired, error_rate=args.error_rate, dispersion_param=args.dispersion)


if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Simulate RNA-seq reads.')

    parser.add_argument('--gtf', metavar='path', type=str,  nargs='+',
                        help='GTF file with gene annotations.  Use with --refernce-fasta.  Alternately, use '
                             '--transcript-fasta.')
    parser.add_argument('--reference-fasta', metavar='path', type=str, nargs='+',
                        help='FASTA files containing reference genome.  Use with --gtf.  Alternately, use '
                             '--transcript-fasta.')
    parser.add_argument('--output-prefix', type=str, default="simulated",
                        help='Prefix to give to all the output file names')
    parser.add_argument('--bio-replicates', type=int, required=False, default=10,
                        help='Number of biological replicates')
    parser.add_argument('--de-fraction', type=float, required=False, default=0.1,
                        help='Fraction of expressed transcript that are differentially expressed between groups')
    parser.add_argument('--error-rate', type=float, required=False, default=0.005,
                        help='Fraction of simulated bases affected by sequencing errors')
    parser.add_argument('--dispersion', type=float, required=False, default=100.0,
                        help='Dispersion parameter affecting amount of biological variation simulated')
    parser.add_argument('--fragment-mean', type=int, required=False, default=250,
                        help='Mean fragment length, drawn from gaussian.  See also: --fragment-sd')
    parser.add_argument('--fragment-sd', type=float, required=False, default=25.0,
                        help='Mean fragment length, drawn from gaussian.  See also: --fragment-sd')
    parser.add_argument('--read-length', type=int, required=False, default=100,  help='Read length')
    parser.add_argument('--paired', action='store_const', const=True, default=False, help='Generate paired-end reads')
    parser.add_argument('--test', action='store_const', const=True, default=False, help='Do unit tests')
    parser.add_argument('--profile', action='store_const', const=True, default=False, help='Print profiling info')
    parser.add_argument('--verbose', action='store_const', const=True, default=False, help='Be talkative')

    args = parser.parse_args()

    # Set up logger
    logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s', datefmt='%m/%d/%y-%H:%M:%S',
                        level=logging.DEBUG if args.verbose else logging.INFO)

    if args.profile:
        import cProfile
        cProfile.run('go(args)')
    else:
        go(args)
