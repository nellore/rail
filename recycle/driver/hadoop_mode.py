"""
hadoop_mode.py

Parameters and objects related to hadoop mode.
"""

def addHadoopModeArgs(parser):
    #
    # Hadoop params
    #
    parser.add_argument(\
        '--hadoop-script', metavar='PATH', type=str, help='Location of Hadoop script')
    parser.add_argument(\
        '--streaming-jar', metavar='PATH', type=str, help='Location of Hadoop streaming jar')
