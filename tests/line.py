#!/usr/bin/env python
"""
line.py

General script for testing consistency of all of Rail's outputs given any
input manifest file and index. Here, consistency is verified using third-party
tools:

  1) bedtools genomecov and bigWigToBedGraph to check that coverage vectors
  from BAM and bigWig are consistent.
  2) samtools mpileup to check that:
     a) mismatch bigWigs and BAMs agree
     b) indel beds and BAMs agree

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
import re

# Global determines whether temp dir should be removed
_should_delete = False

def kill_dir(dir_to_kill):
    """ Removes directory and contents if _should_delete is True

        dir_to_kill: directory to remove

        No return value.
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


def compare_bam_and_bw(sample, working_dir, filter_script, unique=False,
                        bedtools='bedtools', samtools='samtools',
                        bigwigtobedgraph='bigWigToBedGraph'):
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
        filter_script: location of filter.py
        bigwigtobedgraph: path to bigwigtobedgraph executable

        Return value: tuple (True iff no diffs,
                                file containing diffs between bedGraphs)
    """
    print >>sys.stderr, (
                'Comparing {}BAM to bw with bedtools for sample {}...'
            ).format('unique alignments from ' if unique else '', sample)
    alignments_basename = os.path.join(
                            working_dir, 'rail-rna_out',
                            'alignments', 'alignments'
                        )
    coverage_dir = os.path.join(
                            working_dir, 'rail-rna_out',
                            'coverage_bigwigs'
                        )
    bedgraph_from_bam = os.path.join(working_dir, 'from_bam.bedGraph')
    unsorted_bedgraph_from_bw = os.path.join(working_dir,
                                                'from_bw.temp.bedGraph')
    bedgraph_from_bw = os.path.join(working_dir, 'from_bw.bedGraph')
    bam_to_bedgraph_command = ('set -exo pipefail; '
                    '({samtools} view -H {first_bam};'
                    ' for i in {alignments_basename}.{sample}.*.bam;'
                    ' do {samtools} view $i; done) '
                    '| {executable} {filter_script} '
                    '--deletions-to-matches ' +
                    ('--uniques ' if unique else '') +
                    '| {samtools} view -bS - '
                    '| {bedtools} genomecov -ibam - -bga -split '
                    '| sort -k1,1 -k2,2n -k3,3n '
                    '| awk \'$4 != 0\' >{output_bedgraph}'
                ).format(
                    samtools=samtools,
                    first_bam=glob.glob(alignments_basename + '.*.bam')[0],
                    alignments_basename=alignments_basename,
                    filter_script=filter_script,
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
                '{bigwigtobedgraph} {coverage_bigwig} '
                '{unsorted_bedgraph}; '
                'sort -k1,1 -k2,2n -k3,3n {unsorted_bedgraph} '
                '| awk \'$4 != 0\' >{final_bedgraph}; '
                'rm {unsorted_bedgraph};'
            ).format(
                    bigwigtobedgraph=bigwigtobedgraph,
                    coverage_bigwig=os.path.join(
                                coverage_dir,
                                '{sample}{unique}.bw'.format(
                                        sample=sample,
                                        unique=('.unique' if unique else ''),
                                )
                            ),
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
        if e.returncode > 1:
            raise RuntimeError(error_out(e, diff_command))
    if os.path.getsize(diff_output):
        print >>sys.stderr, 'FAIL'
        return (False, diff_output)
    print >>sys.stderr, 'SUCCESS'
    return (True, diff_output)

def compare_bam_and_variants(sample, working_dir, filter_script, unique=False,
                             bedtools='bedtools', samtools='samtools'):
    """ Compares BAMs, bigWigs, and BEDs for a given sample.

        Uses samtools mpileup to obtain mismatches and indels from BAM output
        of Rail-RNA as well as bigWigToBedGraph to obtain bedGraph with 
        mismatches. Gives diff between mpileup indels and indel BEDs from
        Rail-RNA and diff between mpileup mismatches and mismatch bedGraph.

        sample: sample name
        working_dir: directory containing rail-rna_out after running Rail-RNA
        filter_script: location of filter.py
        unique: True iff only unique alignments should be included
        bedtools: path to bedtools executable
        samtools: path to samtools executable

        Return value: tuple (True iff no diffs,
                                file containing diffs between bedGraphs)
    """
    print >>sys.stderr, (
                'Comparing {}BAM to{} mismatch bigWig with '
                'samtools mpileup for sample {}...'
            ).format('unique alignments from ' if unique else '',
                     ' indel BEDs and', sample)
    alignments_basename = os.path.join(
                            working_dir, 'rail-rna_out',
                            'alignments', 'alignments'
                        )
    coverage_dir = os.path.join(
                            working_dir, 'rail-rna_out',
                            'coverage_bigwigs'
                        )
    # Store mismatch tracks from BAM via mpileup
    bam_base = {}
    pileup_command = ('set -exo pipefail; '
                      '({samtools} view -H {first_bam};'
                      ' for i in {alignments_basename}.{sample}.*.bam;'
                      ' do {samtools} view $i; done) '
                      '| {executable} {filter_script} ' + 
                      ('--uniques ' if unique else '') + 
                      '| {samtools} mpileup -B -d 2147483647 '
                      '-Q 0 -q 0 -x --ff -').format(
                        first_bam=glob.glob(alignments_basename + '.*.bam')[0],
                        alignments_basename=alignments_basename,
                        executable=sys.executable,
                        filter_script=filter_script,
                        sample=sample,
                        samtools=samtools
                    )
    pileup_process = subprocess.Popen(pileup_command,
                                        stdout=subprocess.PIPE,
                                        executable='/bin/bash',
                                        shell=True)
    try:
        # For getting bedGraph coordinates right for mismatches
        stretch, last_chrom, last_pos, last_mismatches = (defaultdict(int),
                                                          defaultdict(int),
                                                          defaultdict(int),
                                                          defaultdict(int))
        with open(os.path.join(working_dir, 'from_bam.insertions.bed'), 'w') \
                as bam_insertions, \
              open(os.path.join(working_dir, 'from_bam.deletions.bed'), 'w') \
                    as bam_deletions, \
              open(os.path.join(working_dir, 'from_bam.A.bedGraph'), 'w') \
                    as bam_base['A'], \
              open(os.path.join(working_dir, 'from_bam.C.bedGraph'), 'w') \
                    as bam_base['C'], \
              open(os.path.join(working_dir, 'from_bam.G.bedGraph'), 'w') \
                    as bam_base['G'], \
              open(os.path.join(working_dir, 'from_bam.T.bedGraph'), 'w') \
                    as bam_base['T'], \
              open(os.path.join(working_dir, 'from_bam.N.bedGraph'), 'w') \
                    as bam_base['N']:
            for line in pileup_process.stdout:
                chrom, pos, _, coverage, bases, _ = line.strip().split('\t')
                pos = int(pos) - 1
                pos_str = str(pos)
                insertions = defaultdict(int)
                for insertion in re.finditer(
                                    r'\-[0-9]+([ACGTNacgtn]+)', bases
                                ):
                    insertion_string = insertion.string[
                            insertion.start(1):insertion.end(1)
                        ].upper()
                    insertions[insertion_string] += 1
                for insertion in insertions:
                    print >>bam_insertions, '\t'.join(
                            [chrom, pos_str, pos_str, insertion, coverage]
                        )
                deletions = defaultdict(int)
                for deletion in re.finditer(
                                    r'\-[0-9]+([ACGTNacgtn]+)', bases
                                ):
                    deletion_string = deletion.string[
                            deletion.start(1):deletion.end(1)
                        ].upper()
                    deletions[deletion_string] += 1
                for deletion in deletions:
                    print >>bam_deletions, '\t'.join(
                            [chrom, pos_str, str(pos + len(deletion)),
                                deletion, deletions[deletion]]
                        )
                mismatches = defaultdict(int)
                for mismatch in re.finditer(r'[^0-9]([ACGTNacgtn])'):
                    mismatch_string = mismatch.string[
                            mismatch.start(1):mismatch.end(1)
                        ]
                    mismatches[mismatch_string] += 1
                for mismatch in mismatches:
                    if chrom != last_chrom[mismatch] or (
                            pos > last_pos[mismatch] + 1
                        ) or (
                            mismatches[mismatch] != last_mismatches[mismatch]
                        ):
                        # Write output
                        print >>bam_base[mismatch], '\t'.join(
                                [chrom, last_pos[mismatch] - stretch[mismatch],
                                    last_pos[mismatch],
                                    last_mismatches[mismatch]]
                            )
                        stretch[mismatch] = 0
                    else:
                        stretch[mismatch] += 1
                    last_pos[mismatch] = pos
                    last_mismatches[mismatch] = mismatches[mismatch]
                    last_chrom[mismatch] = chrom
    finally:
        pileup_process.stdout.close()
        pileup_process.wait()
    diff_output = { diff_type : os.path.join(
                        working_dir,
                        '{diff_type}_diffs.{sample}{unique}.txt'.format(
                                diff_type=diff_type,
                                sample=sample,
                                unique=('.unique' if unique else '')
                            )
                    ) for diff_type in ['insertions', 'deletions',
                                        'A', 'C', 'G', 'T', 'N'] }
    for diff_type, from_rail, from_bam in [
                (indel_which,
                 os.path.join(working_dir, 'rail-rna_out',
                                'junctions_and_indels',
                                indel_which + '.' + sample + '.bed'),
                 os.path.join(working_dir, 'from_bam.' + indel_which + '.bed'))
                for indel_which in ['insertions', 'deletions']
            ] + [
            (base,
             os.path.join(working_dir, 'rail-rna_out',
                            'coverage_bigwigs',
                            '{sample}.{base}{unique}.bw'.format(
                                                        sample=sample,
                                                        base=base,
                                                        unique=('.unique'
                                                                if unique
                                                                else '')
                                                    )
                                ),
             os.path.join(working_dir, 'from_bam.' + base + '.bed'))
            for base in 'ACGTN'
        ]:
        diff_command = ('diff <(sort -k1,1 -k2,2n -k3,3n {from_bam}) '
                        '<(cat {from_rail} |{remove_top_line} '
                        'sort -k1,1 -k2,2n -k3,3n) '
                        '>{diffs}').format(
                            bedgraph_from_bam=bedgraph_from_bam,
                            bedgraph_from_bw=bedgraph_from_bw,
                            diffs=os.path.join(working_dir,
                                                diff_type + '.diffs'),
                            remove_top_line=' tail -n +2 |'
                        )
    try:
        subprocess.check_call(
                diff_command, executable='/bin/bash', shell=True
            )
    except subprocess.CalledProcessError as e:
        if e.returncode > 1:
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
    args = parser.parse_args()
    # Get sample names
    with open(args.manifest) as manifest_stream:
        samples = [line.strip().split('\t')[-1] for line in manifest_stream]

    # Perform all work in a temporary directory
    working_dir = tempfile.mkdtemp()
    # Tear down after script is complete
    atexit.register(kill_dir, working_dir)

    rail_dir = os.path.join(
                        os.path.dirname(
                            os.path.dirname(
                                os.path.abspath(__file__)
                                )
                            ),
                        'src'
                    )
    filter_script = os.path.join(
                os.path.dirname(os.path.abspath(__file__)), 'filter.py'
            )
    manifest, bowtie_idx, bowtie2_idx = map(
                                    os.path.abspath,
                                    [args.manifest,
                                        args.bowtie_idx, args.bowtie2_idx]
                                )
    os.chdir(working_dir)
    print >>sys.stderr, 'Running Rail on manifest file "{}"...'.format(
                                                                args.manifest
                                                            )
    try:
        rail_output = subprocess.check_output(
                [sys.executable, rail_dir, 'go', 'local',
                    '-m', manifest, '-x', bowtie_idx, bowtie2_idx,
                    '-d', 'tsv,bed,bam,bw'],
            stderr=sys.stdout)
    except subprocess.CalledProcessError as e:
        error_match = re.search(r'2>(.*\.log)" failed;', e.output)
        if error_match:
            log_file = error_match.string[
                                error_match.start(1):error_match.end(1)
                            ]
            print >>sys.stderr, (
                    'Error occurred running a step. Contents of '
                    'log file {}:'.format(log_file)
                )
            with open(log_file) as log_stream:
                print >>sys.stderr, log_stream.read()
        raise RuntimeError(error_out(e))

    # Check consistency of BAM and bigwig
    for sample, unique in product(samples, (True, False)):
        success, diff_output = compare_bam_and_bw(
                sample, working_dir, filter_script=filter_script,
                unique=unique, bedtools=args.bedtools,
                samtools=args.samtools, bigwigtobedgraph=args.bigwigtobedgraph
            )
        if not success:
            _should_delete = False
            raise RuntimeError(
                    'BAM and bigWig are inconsistent for sample {}. See diffs '
                    'between bedGraphs obtained from either output in '
                    '{}.'
                ).format(sample, diff_output)
