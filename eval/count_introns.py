"""
count_introns.py

Uses bowtie2-inspect to get RNAMEs from transcript fragment index written by
 Rail-RNA and count the number of introns in it. Specification of RNAMEs from
this index is in intron_fasta.py.
"""

if __name__ == '__main__':
    import argparse
    import subprocess

    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--bowtie2-inspect', type=str, required=False,
            default='bowtie2-inspect',
            help='Path to Bowtie 2 executable'
        )
    parser.add_argument('-b', '--basename', type=str, required=True,
            help=('Path to basename of Bowtie 2 index containing transcript '
                  'fragments')
        )

    args = parser.parse_args()

    inspect_process = subprocess.Popen(
                            [args.bowtie2_inspect, '-n', args.basename],
                            stdout=subprocess.PIPE
                        )
    introns = set()
    for line in inspect_process.stdout:
        rname_and_sense, seq_start, subseq_sizes, intron_sizes, _, _ = \
            line.split('\x1d')
        seq_start = int(seq_start)
        subseq_sizes = [int(size) for size in subseq_sizes.split(',')]
        intron_sizes = [int(size) for size in intron_sizes.split(',')]
        assert len(subseq_sizes) == len(intron_sizes) + 1
        start = seq_start
        for i, size in enumerate(subseq_sizes[:-1]):
            introns.add((
                       rname_and_sense,
                       start + size, start + size + intron_sizes[i]
                    )
            )
            start += (size + intron_sizes[i])
    inspect_process.wait()
    print 'intron count: %d' % len(introns)