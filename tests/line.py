#!/usr/bin/env python
"""
line.py

General script for testing consistency of all of Rail's outputs given any
input manifest file and index. Here, consistency is verified using third-party
tools:

  1) bedtools genomecov and bigWigToBedGraph to check that coverage vectors
  from BAM and bigWig are consistent.
  2) https://github.com/lindenb/jvarkit/wiki/SAM2Tsv to check that mismatches
  and indels from BAM and bigWig are consistent.

Requires that bedtools, samtools, and SAM2Tsv are installed. Performs all tests
in local mode for now.
"""
import tempfile
import atexit
import shutil
import os
import subprocess
import traceback
import sys
import glob
from itertools import product
from shutil import rmtree

# Global determines whether temp dir should be removed
_should_delete = True

def kill_dir(dir_to_kill):
    """ Removes directory and contents if _should_delete is True

        dir_to_kill: directory to remove
    """
    if _should_delete:
        rmtree(dir_to_kill)

def error_out(e, command=None):
    """ Returns message for last exception

        e: Exception object

        Return value: string reporting error
    """
    return 'Error "{}" encountered{}.\n{}'.format(
                                                    e.message,
                                                    ' executing {}.'.format(
                                                                        command
                                                                    )
                                                    if command is not None
                                                    else '',
                                                    traceback.format_exc()
                                                )


def compare_bam_and_bw(sample, working_dir, unique=False,
                        bedtools='bedtools', samtools='samtools',
                        bedgraphtobigwig='bedGraphToBigWig'):
    """ Compares BAM and bigwig for a given sample.

        Uses bedtools genomecov to obtain bedGraph encoding coverage vector
        from BAM output of Rail-RNA as well as bigWigToBedGraph to obtain
        bedGraph encoding coverage vector from bed output of Rail-RNA.
        Gives diff between bedGraphs.

        sample: sample name
        working_dir: directory containing rail-rna_out after running Rail-RNA
        unique: True iff only unique alignments should be included
        bedtools: path to bedtools executable
        samtools: path to samtools executable
        bedgraphtobigwig: path to bedgraphtobigwig executable

        Return value: tuple (True iff no diffs,
                                file containing diffs between bedGraphs)
    """
    print >>sys.stderr, (
                'Comparing {}BAM to bw with bedtools for sample {}...'
            ).format('unique alignments from ' if unique else '')
    alignments_basename = os.path.join(working_dir, 'alignments')
    bedgraph_from_bam = os.path.join(working_dir, 'from_bam.bedGraph')
    unsorted_bedgraph_from_bw = os.path.join(working_dir,
                                                'from_bw.temp.bedGraph')
    bedgraph_from_bw = os.path.join(working_dir, 'from_bw.bedGraph')
    bam_to_bedgraph_command = (
                    '{samtools} view -H {first_bam}; '
                    'for i in {alignments_basename}.{sample}.bam; '
                    'do {samtools} view $i; done '
                    '| {executable} filter.py --deletions-to-matches ' +
                    ('--uniques ' if unique else '') +
                    '| {samtools} view -bS - '
                    '| {bedtools} genomecov -ibam - -bga -split '
                    '| sort -k1,1 -k2,2n >{output_bedgraph}'
                ).format(
                    samtools=samtools,
                    first_bam=glob.glob(os.path.join(
                                                working_dir,
                                                'alignments.*.bam'
                                            ))[0],
                    alignments_basename=os.path.join(
                                                working_dir,
                                                'alignments'
                                            ),
                    sample=sample,
                    executable=sys.executable,
                    bedtools=bedtools,
                    output_bedgraph=bedgraph_from_bam
            )
    try:
        subprocess.check_output(
                bam_to_bedgraph_command, executable='/bin/bash', shell=True
            )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(error_out(e, bam_to_bedgraph_command))
    bw_to_bedgraph_command = (
                'set -exo pipefail; '
                '{bigwigtobedgraph} {sample}{unique}.bw '
                '{unsorted_bedgraph}; '
                'sort -k1,1 -k2,2n {unsorted_bedgraph} >{final_bedgraph}; '
                'rm {unsorted_bedgraph};'
            ).format(
                    bigwigtobedgraph=bigwigtobedgraph,
                    sample=sample,
                    unique=('.unique' if unique else ''),
                    unsorted_bedgraph=unsorted_bedgraph_from_bw,
                    final_bedgraph=bedgraph_from_bw
                )
    try:
        subprocess.check_call(
                bw_to_bedgraph_command, executable='/bin/bash', shell=True
            )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(error_out(e, bw_to_bedgraph_command))
    diff_output = os.path.join(
                        working_dir,
                        'coverage_diffs.{sample}{unique}.txt'.format(
                                sample=sample,
                                unique=('.unique' if unique else '')
                            )
                    )
    diff_command = ('diff {bedgraph_from_bam} '
                    '{bedgraph_from_bw} >{diffs}').format(
                        bedgraph_from_bam=bedgraph_from_bam,
                        bedgraph_from_bw=bedgraph_from_bw,
                        diffs=diff_output
                    )
    try:
        subprocess.check_call(
                diff_command, executable='/bin/bash', shell=True
            )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(error_out(e, diff_command))
    if os.path.getsize(diff_output):
        print >>sys.stderr, 'FAIL'
        return (False, diff_output)
    print >>sys.stderr, 'SUCCESS'
    return (True, diff_output)

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(
                description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter
            )
    # Add command-line arguments
    parser.add_argument('--manifest', '-m', type=str, required=True,
            help='Rail-RNA test manifest file'
        )
    parser.add_argument('--bowtie-idx', type=str, required=True,
            help='Bowtie index basename'
        )
    parser.add_argument('--bowtie2-idx', type=str, required=True,
            help='Bowtie2 index basename'
        )
    parser.add_argument('--bedtools', type=str, required=False,
            default='bedtools',
            help='Bedtools executable; known to work on v'
        )
    parser.add_argument('--samtools', type=str, required=False,
            default='samtools',
            help='SAMTools executable; known to work on v'
        )
    parser.add_argument('--bigwigtobedgraph', type=str, required=False,
            default='bigWigToBedGraph',
            help='bigWigToBedGraph executable; unversioned Kent Tool'
        )
    parser.add_argument('--sam2tsv', type=str, required=False,
            default='SAM2Tsv.jar',
            help='SAM2Tsv jar; known to work on v'
        )
    args = parser.parse_args()
    # Get sample names
    with open(args.manifest) as manifest_stream:
        samples = [line.strip().split('\t')[-1] for line in manifest_stream]

    # Perform all work in a temporary directory
    working_dir = tempfile.mkdtemp()
    # Tear down after script is complete
    atexit.register(kill_dir, working_dir)

    rail_dir = os.path.join(
                        os.path.dirname(os.path.dirname(__file__)), 'src'
                    )
    os.chdir(working_dir)
    print >>sys.stderr, 'Running Rail on manifest file "{}"...'.format(
                                                                args.manifest
                                                            )
    try:
        subprocess.check_call([sys.executable, rail_dir, 'go', 'local',
                                 '-m', args.manifest, '-x', bowtie_idx,
                                 bowtie2_idx, '-d',
                                 'tsv,bed,bam,bw'], stderr=sys.stderr,
                                 shell=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(error_out(e))

    # Check consistency of BAM and bigwig
    for sample, unique in product(samples, (True, False)):
        success, diff_output = compare_bam_and_bw(
                sample, working_dir, unique=unique, bedtools=args.bedtools,
                samtools=args.samtools, bedgraphtobigwig=args.bedgraphtobigwig
            )
        if not success:
            raise RuntimeError(
                    'BAM and bigWig are inconsistent for sample {}. See diffs '
                    'between bedGraphs obtained from either output in '
                    '{}.'
                ).format(sample, diff_output)
