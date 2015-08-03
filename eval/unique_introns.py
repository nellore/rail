"""
unique_introns.py

Measures precision and recall of introns that are unique to each of the 20
samples selected at random (in create_single_sample_sim_commands.py) from
112 simulated GEUVADIS-based samples. Assumes a directory structure for output
created by either grab_112_sim_results.sh or 

A file Aligned.out.sam is assumed to be in output directories for all aligners
but Rail-RNA and TopHat 2. TopHat 2's output file is assumed to be
accepted_hits.bam, and Rail's output is assumed to span all the bam files in
its output directory.
"""

import os
from collections import defaultdict
import sys

def introns_from_bed(bed, index):
    """ Converts BED to dictionary that maps RNAMES to sets of introns.

        bed: input BED file characterizing splice junctions
        index: index of bed file; used for tracking which bed is which
            when using multiple cores

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
    return (index, introns)

def introns_from_sam_stream(sam_stream):
    """ Writes output that maps QNAMES to exon-exon junctions overlapped.

        sam_stream: where to find retrieved alignments in SAM form

        Return value: a dictionary. Each key is an RNAME, typically a
            chromosome, and its corresponding value is a set of tuples, each
            denoting an intron on RNAME. Each tuple is of the form
            (start position, end position).
    """
    introns = defaultdict(set)
    for line in sam_stream:
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
            #md = [token[5:] for token in tokens if token[:5] == 'MD:Z:'][0]
            _, _, introns_to_add, _ = indels_introns_and_exons(cigar,
                                        dummy_md_index(cigar), pos, seq)
            for intron in introns_to_add:
                introns[rname].add(intron[:2])
                key = (rname,) + intron[:2]
        except IndexError:
            print >>sys.stderr, ('Error found on line: ' + line)
            raise
    return introns

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
            required=False,
            default=None,
            help='Path to directory containing Flux Simulator BEDs'
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
    # Use multiple cores
    pool = multiprocessing.Pool(multiprocessing.cpu_count() - 1)
    for i, bed_path in enumerate(bed_paths):
        with open(sample_sim_path) as bed_stream:
            for intron in introns_from_bed(bed_path):
                true_introns[intron].append(i)
    true_unique_introns = defaultdict(set)
    for intron in true_introns:
        if len(true_introns[intron]) == 1:
            true_unique_introns[true_introns[intron][0]].add(intron)
    if 'star' in args.aligners:
        for mode in ['ann_paired_1pass', 'ann_paired_2pass',
                     'noann_paired_1pass', 'nogen_noann_paired_2pass']:

    if args.true_introns_bed_dir is not None:
        # Read Flux BEDs
        true_introns = set()
        introns = set([(strand[:-1], pos, end_pos) for (strand, pos, end_pos) in introns])
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