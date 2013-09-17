"""
config.py

Parameters that are independent of mode (local / hadoop / emr) and independent
of the particular tool being run.
"""

import os

def addConfigArgs(parser):
    #
    # Modes of operation
    #
    parser.add_argument(\
        '--local', action='store_const', const=True, help='Run in local mode.')
    parser.add_argument(\
        '--hadoop', action='store_const', const=True, help='Run in Hadoop mode.')
    parser.add_argument(\
        '--emr', action='store_const', const=True, help='Run in Elastic MapReduce mode.')
    parser.add_argument(\
        '--test', action='store_const', const=True, help='Just test to see if requisite files are around.')
    
    #
    # Basic plumbing
    #
    parser.add_argument(\
        '--input', metavar='PATH', type=str, required=False, help='URL for input directory')
    parser.add_argument(\
        '--output', metavar='PATH', type=str, required=True, help='URL for output directory')
    parser.add_argument(\
        '--reference', metavar='PATH', type=str, required=False, help='URL for reference archive')
    parser.add_argument(\
        '--igenomes', action='store_const', const=True, help='--reference specifies an iGenomes archive/directory, not a Myrna 2 reference')
    parser.add_argument(\
        '--intermediate', metavar='PATH', type=str, help='URL for intermediate files')
    parser.add_argument(\
        '--keep-intermediates', action='store_const', const=True, help='Keep intermediate files in addition to final output files.')
    parser.add_argument(\
        '--dry-run', action='store_const', const=True, help='Just generate script for launching EMR cluster, but don\'t launch it.')
    
    #
    # Preprocessing params
    #
    parser.add_argument(\
        '--preprocess-output', metavar='PATH', type=str, help='Put output from preprocessing step here')
    parser.add_argument(\
        '--preprocess-compress', metavar='gzip|none|bzip2', type=str, default='gzip', help='Type of compression to use for preprocessing output.')
    
    #
    # Other params
    #
    parser.add_argument(\
        '--verbose', action='store_const', const=True, help='Print lots of info to stderr.')
    parser.add_argument(\
        '--version', action='store_const', const=True, help='Just print version information and quit.')
    parser.add_argument(\
        '--set-version', metavar='VER', type=str, help='Force Tornado to use a particular version.')

class GenericConfig(object):
    def __init__(self, args, out):
        self.preprocCompress = args.preprocess_compress
        self.out = out.toUpperUrl()
        if out.isLocal():
            self.out = os.path.abspath(self.out)
