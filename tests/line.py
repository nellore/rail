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

Requires that bedtools and SAM2Tsv are both installed. Performs tests locally
for now.
"""
import tempfile
import atexit
import shutil
import os
import subprocess

if __name__ == '__main__':
	import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
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
    parser.add_argument('--sam2tsv', type=str, required=False,
            default='SAM2Tsv.jar',
            help='SAM2Tsv jar; known to work on v'
        )
    # Get sample names
    with open(args.manifest) as manifest_stream:
        sample_names = [line.strip().split('\t')[-1]
                            for line in manifest_stream]

    # Perform all work in a temporary directory
    working_dir = tempfile.mkdtemp()
    # Tear down after script is complete
    atexit.register(shutil.rmtree, working_dir)

    rail_dir = os.path.join(
                        os.path.dirname(os.path.dirname(__file__)), 'src'
                    )
    os.chdir(working_dir)
    print >>sys.stderr, 'Running Rail on manifest file "{}"...'.format(
                                                                args.manifest
                                                            )
    subprocess.check_output(['python', rail_dir, 'go', 'local',
                                '-m', args.manifest, '-x', bowtie_idx,
                                bowtie2_idx])
    