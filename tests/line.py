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
  3) samtools view to verify that sequences from BAM are the same as
    intermediate sequences obtained after preprocessing

Requires that bedtools an samtools are installed. Performs all tests in local
mode for now.

TODO: 1) Check that junction BEDs agree with cross-sample junction outputs
    from _both_ first and second-pass alignment
    2) Check that indel BEDs agree with cross-sample indel TSVs
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
from collections import defaultdict
import string
import gzip

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

def diff_exit(output_file):
    """ Returns file contents if using Travis CI; else prints filename.

        output_file: path to output file

        Return value: string with file contents if using Travis; else filename
    """
    if 'HAS_JOSH_K_SEAL_OF_APPROVAL' in os.environ:
        # Inside joke? that means we're in Travis CI environment
        with open(output_file) as output_stream:
            return '\n'.join(
                    ['Diffs from {}:'.format(output_file),
                     output_stream.read()]
                )
    return 'See diffs in {}.'.format(output_file)


def compare_bam_and_raw_data(sample, working_dir, samtools='samtools'):
    """ Tests whether read sequences in raw data and BAMs are the same.

        Accounts for how FLAG in BAM may indicate that an alignment record's
        SEQ field is reverse-complemented. Here, raw data is the output of 
        preprocess.py.

        sample: sample name
        working_dir: directory containing rail-rna_out after running Rail-RNA
        samtools: path to samtools executable

        Return value: tuple (True iff no diffs,
                                file containing diffs between name-sequence
                                combinations)
    """
    print >>sys.stderr, (
            'Comparing preprocessed SEQs with BAM SEQs for sample {}...'
        ).format(sample)
    # Output files before sorting and diffing
    from_prep = os.path.join(working_dir, 'from_prep.tsv')
    from_bam = os.path.join(working_dir, 'from_bam.tsv')
    # For reverse complementing
    translation_table = string.maketrans('ATCGN', 'TAGCN')
    with open(from_prep, 'w') as output_stream:
        for prep_file in glob.glob(os.path.join(
                                working_dir, 'rail-rna_logs',
                                'preprocess', 'push', '*.gz'
                            )):
            with gzip.open(prep_file) as prep_stream:
                for line in prep_stream:
                    (seq, reverse_complemented,
                            qname, _) = line.strip().split('\t')
                    qname, _, current_sample = qname.split('\x1d')
                    if current_sample != sample: continue
                    if reverse_complemented == '1':
                        seq = seq[::-1].translate(translation_table)
                    if qname.endswith('/1') or qname.endswith('/2'):
                        qname = qname[:-2]
                    print >>output_stream, '\t'.join([qname, seq])
    bams = os.path.join(
                        working_dir, 'rail-rna_out', 'alignments',
                        'alignments.' + sample + '.*.bam'
                    )
    samtools_command = ('set -exo pipefail; '
                        'for i in {bams}; '
                        'do {samtools} view $i; done').format(
                            bams=bams,
                            samtools=samtools
                        )
    samtools_process = subprocess.Popen(samtools_command,
                                        stdout=subprocess.PIPE,
                                        executable='/bin/bash',
                                        shell=True)
    try:
        with open(from_bam, 'w') as output_stream:
            for line in samtools_process.stdout:
                (qname, flag, _, _, _, _,
                    _, _, _, seq) = line.strip().split('\t')[:10]
                if int(flag) & 16:
                    seq = seq[::-1].translate(translation_table)
                if qname.endswith('/1') or qname.endswith('/2'):
                    qname = qname[:-2]
                print >>output_stream, '\t'.join([qname, seq])
    finally:
        samtools_process.stdout.close()
        samtools_process.wait()
    diff_output = os.path.join(working_dir, 'seq_diffs.' + sample + '.txt')
    diff_command = ('diff <(sort {from_bam}) '
                    '<(sort {from_prep}) >{diffs}').format(
                        from_bam=from_bam,
                        from_prep=from_prep,
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

def compare_bam_and_bw(sample, working_dir, filter_script, unique=False,
                        bedtools='bedtools', samtools='samtools',
                        bigwigtobedgraph='bigWigToBedGraph'):
    """ Compares BAM and bigwig for a given sample.

        Uses bedtools genomecov to obtain bedgraph encoding coverage vector
        from BAM output of Rail-RNA as well as bigWigToBedGraph to obtain
        bedgraph encoding coverage vector from bed output of Rail-RNA.
        Gives diff between bedgraphs.

        sample: sample name
        working_dir: directory containing rail-rna_out after running Rail-RNA
        unique: True iff only unique alignments should be included
        bedtools: path to bedtools executable
        samtools: path to samtools executable
        filter_script: location of filter.py
        bigwigtobedgraph: path to bigwigtobedgraph executable

        Return value: tuple (True iff no diffs,
                                file containing diffs between bedgraphs)
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
    bedgraph_from_bam = os.path.join(working_dir, 'from_bam.bedgraph')
    unsorted_bedgraph_from_bw = os.path.join(working_dir,
                                                'from_bw.temp.bedgraph')
    bedgraph_from_bw = os.path.join(working_dir, 'from_bw.bedgraph')
    '''awk below eliminates all coverage-0 lines because bedGraphToBigWig
    will leave out empty chromosomes, making for diffs on comparison.'''
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

def compare_bam_and_variants(sample, working_dir, filter_script, genome,
                             unique=False, bedtools='bedtools',
                             samtools='samtools',
                             bigwigtobedgraph='bigWigToBedGraph'):
    """ Compares BAMs, bigWigs, and BEDs for a given sample.

        Uses samtools mpileup to obtain mismatches and indels from BAM output
        of Rail-RNA as well as bigWigToBedGraph to obtain bedgraph with 
        mismatches. Gives diff between mpileup indels and indel BEDs from
        Rail-RNA and diff between mpileup mismatches and mismatch bedgraph.
        Note that indels aren't compared when unique=True because there are 
        no BEDs recording indels from uniquely aligning reads in Rail's output.

        sample: sample name
        working_dir: directory containing rail-rna_out after running Rail-RNA
        filter_script: location of filter.py
        genome: path to faidx-indexed reference fasta
        unique: True iff only unique alignments should be included
        bedtools: path to bedtools executable
        samtools: path to samtools executable

        Return value: tuple (True iff no diffs,
                                file containing diffs between bedgraphs)
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
    # Get chromosome order
    first_bam_view = '{samtools} view -H {first_bam}'.format(
                        samtools=samtools,
                        first_bam=glob.glob(alignments_basename + '.*.bam')[0]
                    )
    chroms = subprocess.check_output(
                ('{first_bam_view} '
                 '| grep "@SQ" | cut -f2 | cut -d":" -f2').format(
                       first_bam_view=first_bam_view
                    ), executable='/bin/bash', shell=True
            )
    bams = ['{alignments_basename}.{sample}.{chrom}.bam'.format(
                                    alignments_basename=alignments_basename,
                                    sample=sample,
                                    chrom=chrom
                                ) for chrom in chroms.strip().split('\n')]
    # Eliminate chroms to which reads did not align
    bams = ' '.join([bam for bam in bams if os.path.exists(bam)])
    # Store mismatch tracks from BAM via mpileup
    bam_base = {}
    pileup_command = ('set -exo pipefail; '
                      '({first_bam_view};'
                      ' for i in {bams};'
                      ' do {samtools} view $i; done) '
                      '| {executable} {filter_script} ' + 
                      ('--uniques ' if unique else '') + 
                      '| {samtools} mpileup -B -d 2147483647 -A -B -o 0 -F 0 '
                      '-Q 0 -q 0 -x -f {genome} ' +
                      ('-I ' if unique else '') + 
                      '--ff UNMAP,SECONDARY -').format(
                        first_bam_view=first_bam_view,
                        bams=bams,
                        executable=sys.executable,
                        filter_script=filter_script,
                        samtools=samtools,
                        genome=genome
                    )
    pileup_process = subprocess.Popen(pileup_command,
                                        stdout=subprocess.PIPE,
                                        executable='/bin/bash',
                                        shell=True)
    try:
        # For getting bedgraph coordinates right for mismatches
        stretch, last_chrom, last_pos, last_mismatches = (defaultdict(int),
                                                          defaultdict(int),
                                                          defaultdict(int),
                                                          defaultdict(int))
        with open(os.path.join(working_dir, 'from_bam.insertions.bed'), 'w') \
                as bam_insertions, \
              open(os.path.join(working_dir, 'from_bam.deletions.bed'), 'w') \
                    as bam_deletions, \
              open(os.path.join(working_dir, 'from_bam.A.bedgraph'), 'w') \
                    as bam_base['A'], \
              open(os.path.join(working_dir, 'from_bam.C.bedgraph'), 'w') \
                    as bam_base['C'], \
              open(os.path.join(working_dir, 'from_bam.G.bedgraph'), 'w') \
                    as bam_base['G'], \
              open(os.path.join(working_dir, 'from_bam.T.bedgraph'), 'w') \
                    as bam_base['T'], \
              open(os.path.join(working_dir, 'from_bam.N.bedgraph'), 'w') \
                    as bam_base['N']:
            for line in pileup_process.stdout:
                chrom, pos, _, _, bases, _ = line.strip().split('\t')
                pos = int(pos) - 1
                pos_str = str(pos)
                mismatches = []
                if not unique:
                    insertions = defaultdict(int)
                    for insertion in re.finditer(
                                        r'\+([0-9]+)([ACGTNacgtn]+)', bases
                                    ):
                        insertion_string = insertion.string[
                                insertion.start(2):insertion.start(2)
                                + int(insertion.string[
                                            insertion.start(1):insertion.end(1)
                                        ])
                            ].upper()
                        mismatches.append(
                                    insertion.string[
                                        insertion.start(2):insertion.end(2)
                                    ][len(insertion_string):]
                                )
                        insertions[insertion_string] += 1
                    for insertion in insertions:
                        print >>bam_insertions, '\t'.join(map(str,
                                    [chrom, pos_str, pos_str, insertion,
                                        insertions[insertion_string]]
                                )
                            )
                    deletions = defaultdict(int)
                    for deletion in re.finditer(
                                        r'\-([0-9]+)([ACGTNacgtn]+)', bases
                                    ):
                        deletion_string = deletion.string[
                                deletion.start(2):deletion.start(2)
                                + int(deletion.string[
                                            deletion.start(1):deletion.end(1)
                                        ])
                            ].upper()
                        mismatches.append(
                                    deletion.string[
                                        deletion.start(2):deletion.end(2)
                                    ][len(deletion_string):]
                                )
                        deletions[deletion_string] += 1
                    for deletion in deletions:
                        print >>bam_deletions, '\t'.join(map(str, 
                                    [chrom, pos + 1, pos + 1 + len(deletion),
                                        deletion, deletions[deletion]]
                                )
                            )
                '''Residual mismatches are preceded by neither deletion
                nor insertion'''
                for mismatch in re.finditer(
                        r'(^|,|\.|\^.|<|>|\$)([ACGTNacgtn]+)', bases
                    ):
                    mismatches.append(mismatch.string[
                            mismatch.start(2):mismatch.end(2)
                        ].upper())
                mismatches = ''.join(mismatches)
                mismatches = { base : mismatches.count(base)
                                for base in 'ATCGN' }
                for mismatch in mismatches:
                    if not mismatches[mismatch]: continue
                    if stretch[mismatch] and (
                            chrom != last_chrom[mismatch] or (
                                pos > last_pos[mismatch] + 1
                            ) or (
                                mismatches[mismatch]
                                    != last_mismatches[mismatch]
                            )
                        ):
                        # Write output
                        print >>bam_base[mismatch], '\t'.join(map(str,
                                    [last_chrom[mismatch],
                                        last_pos[mismatch] + 1
                                            - stretch[mismatch],
                                        last_pos[mismatch] + 1,
                                        last_mismatches[mismatch]]
                                )
                            )
                        stretch[mismatch] = 1
                    else:
                        stretch[mismatch] += 1
                    last_mismatches[mismatch] = mismatches[mismatch]
                    last_pos[mismatch] = pos
                    last_chrom[mismatch] = chrom
            # Last stretch has yet to be printed
            for mismatch in last_mismatches:
                print >>bam_base[mismatch], '\t'.join(map(str,
                                    [last_chrom[mismatch],
                                        last_pos[mismatch] + 1
                                            - stretch[mismatch],
                                        last_pos[mismatch] + 1,
                                        last_mismatches[mismatch]]
                                )
                            )
    finally:
        pileup_process.stdout.close()
        pileup_process.wait()
    # Store mismatch tracks from bw via bigWigToBedGraph
    unsorted_bedgraph_from_bw = os.path.join(working_dir,
                                             'from_bw.temp.bedgraph')
    for base in 'ATCGN':
        bedgraph_from_bw = os.path.join(working_dir,
                                        'from_bw.{base}.bedgraph'.format(
                                            base=base
                                        ))
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
                                    '{sample}.{base}{unique}.bw'.format(
                                        sample=sample,
                                        unique=('.unique' if unique else ''),
                                        base=base
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
    diff_output = { diff_type : os.path.join(
                        working_dir,
                        '{diff_type}_diffs.{sample}{unique}.txt'.format(
                                diff_type=diff_type,
                                sample=sample,
                                unique=('.unique' if unique else '')
                            )
                    ) for diff_type in ['insertions', 'deletions',
                                        'A', 'C', 'G', 'T', 'N'] }
    for diff_type, from_rail, from_bam in ([
                (indel_which,
                 os.path.join(working_dir, 'rail-rna_out',
                                'junctions_and_indels',
                                indel_which + '.' + sample + '.bed'),
                 os.path.join(working_dir, 'from_bam.' + indel_which + '.bed'))
                for indel_which in ['insertions', 'deletions']
            ] if not unique else []) + [
            (base,
             os.path.join(working_dir, 'from_bw.' + base + '.bedgraph'),
             os.path.join(working_dir, 'from_bam.' + base + '.bedgraph'))
            for base in 'ACGTN'
        ]:
        print >>sys.stderr, 'Comparing {}...'.format(diff_type)
        diff_output = os.path.join(
                        working_dir,
                        '{diff_type}_diffs.{sample}{unique}.txt'.format(
                                diff_type=diff_type,
                                sample=sample,
                                unique=('.unique' if unique else '')
                            )
                    )
        diff_command = ('set -exo pipefail; '
                        'diff <(sort -k1,1 -k2,2n -k3,3n {from_bam}) '
                        '<(cat {from_rail} |{remove_top_line} '
                        'sort -k1,1 -k2,2n -k3,3n) '
                        '>{diffs}').format(
                            from_bam=from_bam,
                            from_rail=from_rail,
                            diffs=diff_output,
                            # Track line needs removal for indel beds
                            remove_top_line=(' tail -n +2 |'
                                    if diff_type in ['insertions', 'deletions']
                                    else '')
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
    # Don't return a diff filename; no need
    return (True, '')

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
    parser.add_argument('--genome', type=str, required=True,
            help='faidx-indexed reference fasta'
        )
    parser.add_argument('--bedtools', type=str, required=False,
            default='bedtools',
            help='Bedtools executable; known to work on v2.25'
        )
    parser.add_argument('--samtools', type=str, required=False,
            default='samtools',
            help='SAMTools executable; known to work on v1.3'
        )
    parser.add_argument('--bigwigtobedgraph', type=str, required=False,
            default='bigWigToBedGraph',
            help='bigWigToBedGraph executable; unversioned Kent Tool'
        )
    args = parser.parse_args()
    # Check that genome exists because SAMTools doesn't do this!
    if not os.path.isfile(args.genome):
        raise IOError(
                'Reference FASTA "{}" does not exist.'.format(args.genome)
            )
    # Get absolute path to FASTA since we'll be changing dir
    args.genome = os.path.abspath(args.genome)
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
                [sys.executable, rail_dir, 'go', 'local', '-p', '1',
                    '-m', manifest, '-x', bowtie_idx, bowtie2_idx,
                    '-d', 'tsv,bed,bam,bw', '--gzip-intermediates',
                    '--keep-intermediates'],
                stderr=subprocess.STDOUT
            )
    except subprocess.CalledProcessError as e:
        error_match = re.search(r'2> >\(tee (.*\.log)', e.output)
        if error_match:
            log_file = error_match.string[
                                error_match.start(1):error_match.end(1)
                            ]
            print >>sys.stderr, (
                    'Error occurred running a step. Contents of '
                    'log file {}, where error was found:'.format(log_file)
                )
            with open(log_file) as log_stream:
                print >>sys.stderr, log_stream.read()
            for other_log_file in glob.glob(
                        os.path.join(os.path.dirname(log_file), '*.log')
                    ):
                if other_log_file != log_file:
                    print >>sys.stderr, 'Contents of log file {}:'.format(
                                                            other_log_file
                                                        )
                    with open(other_log_file) as log_stream:
                        print >>sys.stderr, log_stream.read()

        raise RuntimeError(error_out(e))

    for sample, unique in product(samples, (False, True)):
        if not unique:
            # Compare seqs between BAM and preprocessed intermediates
            success, diff_output = compare_bam_and_raw_data(
                    sample, working_dir, samtools=args.samtools
                )
            if not success:
                _should_delete = False
                raise RuntimeError(
                        ('Sequences from BAM and preprocessed inputs do not '
                         'agree for sample {}. {}').format(
                                sample, diff_exit(diff_output)
                            )
                    )
        # Check consistency of BAM and bigwig
        success, diff_output = compare_bam_and_bw(
                sample, working_dir, filter_script=filter_script,
                unique=unique, bedtools=args.bedtools,
                samtools=args.samtools, bigwigtobedgraph=args.bigwigtobedgraph
            )
        if not success:
            _should_delete = False
            raise RuntimeError(
                    ('BAM and bigWig are inconsistent for sample {}. {}'
                     ).format(sample, diff_exit(diff_output))
                )
        # Check consistency of BAM and variant files
        success, diff_output = compare_bam_and_variants(
                sample, working_dir, filter_script=filter_script,
                genome=args.genome, unique=unique, bedtools=args.bedtools,
                samtools=args.samtools, bigwigtobedgraph=args.bigwigtobedgraph,
            )
        if not success:
            _should_delete = False
            raise RuntimeError(
                    'BAM and variant BEDs/bigWigs are inconsistent for sample '
                    '{}. {}'.format(sample, diff_exit(diff_output))
                )

    print >>sys.stderr, 'All line tests finished successfully.'
