#!/usr/bin/env python
"""
bowtie.py
Part of Rail-RNA

Contains Bowtie-related command-line parameters common to steps. Also has
a workaround for a Bowtie 2 bug.
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
