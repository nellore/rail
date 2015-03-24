"""
local_mode.py

Parameters and objects related to local mode.
"""

def addLocalModeArgs(parser):
    parser.add_argument(\
        '--num-processes', metavar='INT', type=int, default=1, help='Max # simultaneous mappers/reducers to run.')
    parser.add_argument(\
        '--dont-overwrite', action='store_const', const=True, help='Abort rather than overwrite any existing directories.')
    parser.add_argument(\
        '--keep-all', action='store_const', const=True, help='Keep all intermediate directories.')

class LocalConfig(object):
    def __init__(self, args):
        self.numProcesses = args.num_processes
        self.force = not args.dont_overwrite
        self.keepAll = args.keep_all
