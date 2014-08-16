#!/usr/bin/env python
"""
bowtie.py
Part of Rail-RNA

Contains Bowtie-related command-line parameters common to steps. Also has
a function for parsing certain Bowtie2 arguments.
"""

def add_args(parser):
    parser.add_argument(\
        '--bowtie-exe', metavar='EXE', type=str, required=False,
        default='bowtie',
        help='Path to executable for Bowtie')
    parser.add_argument(\
        '--bowtie-build-exe', metavar='EXE', type=str, required=False,
        default='bowtie-build',
        help='Path to executable for Bowtie-build')
    parser.add_argument(\
        '--bowtie-idx', metavar='INDEX', type=str, required=False,
        default='',
        help='Path to Bowtie index. Specify its basename.')
    parser.add_argument(\
        '--bowtie2-exe', metavar='EXE', type=str, required=False,
        default='bowtie2',
        help='Path to executable for Bowtie2')
    parser.add_argument(\
        '--bowtie2-build-exe', metavar='EXE', type=str, required=False,
        default='bowtie2-build',
        help='Path to executable for Bowtie2-build')
    parser.add_argument(\
        '--bowtie2-idx', metavar='INDEX', type=str, required=False,
        default='',
        help='Path to Bowtie2 index. Specify its basename.')

def parsed_bowtie_args(bowtie2_args):
    """ Parses Bowtie2 args and returns relevant information about reporting.

        bowtie2_args: string of Bowtie 2 arguments to parse. Parsed arguments
            include -k and -a.

        Return value: tuple (-1 iff all alignments are to be reported; else the
            number of alignments to be reported, --seed parameter,
            --non-deterministic parameter)
    """
    import argparse
    bowtie_parser = argparse.ArgumentParser()
    '''By default, report primary alignment; this is regarded as '-k 1'. Note
    that Bowtie2 does not guarantee that the best alignment is reported when
    -k 1 is invoked, but Rail does here since it has all alignments to work
    with.'''
    bowtie_parser.add_argument('-k', type=int, required=False, default=1)
    bowtie_parser.add_argument('-a', action='store_const', const=True,
            default=False
        )
    bowtie_parser.add_argument('--seed', type=int, required=False, default=0)
    bowtie_parser.add_argument('--non-deterministic', action='store_const',
            const=True, default=False
        )
    if bowtie2_args is None: bowtie2_args = ''
    split_args = bowtie2_args.split(' ')
    parsed_args = bowtie_parser.parse_known_args(split_args)[0]
    try:
        # If both -k and -a are present, last argument takes precedence
        if split_args.index('-a') > split_args.index('-k'):
            return -1, parsed_args.seed, parsed_args.non_deterministic
        else:
            return (
                    parsed_args.k, parsed_args.seed,
                    parsed_args.non_deterministic
                )
    except ValueError:
        # Both -a and -k are not present
        pass
    if parsed_args.a:
        return -1, parsed_args.seed, parsed_args.non_deterministic
    return parsed_args.k, parsed_args.seed, parsed_args.non_deterministic