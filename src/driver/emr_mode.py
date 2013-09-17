"""
emr_mode.py

Parameters and objects related to emr mode.
"""

def addEmrModeArgs(parser):
    #
    # Elastic MapReduce params
    #
    parser.add_argument(\
        '--emr-script', metavar='PATH', type=str, help='Path to Amazon elastic-mapreduce script')
    parser.add_argument(\
        '--emr-local-dir', metavar='PATH', type=str, help='Path to a local directory on the EMR nodes where the reference archive and Tornado scripts will be copied.')
    parser.add_argument(\
        '--hadoop-version', metavar='VERS', type=str, help='Hadoop version number to use')
    parser.add_argument(\
        '--credentials', metavar='PATH', type=str, help='Amazon Elastic MapReduce credentials file')
    parser.add_argument(\
        '--alive', action='store_const', const=True, help='Keep cluster alive after it finishes job.')
    parser.add_argument(\
        '--no-memory-intensive', action='store_const', const=True, help='Do not set memory-intensive mode on the cluster.')
    parser.add_argument(\
        '--no-add-swap', action='store_const', const=True, help='Do not add swap memory to the cluster.')
    parser.add_argument(\
        '--enable-speculative', action='store_const', const=True, help='Enable Hadoop speculative execution for EMR runs (not recommended).')
    parser.add_argument(\
        '--name', metavar='STR', type=str, help='Amazon Elastic MapReduce job name.')
    parser.add_argument(\
        '--no-emr-debug', action='store_const', const=True, help='Don\'t enable EMR debugging functions.')
    parser.add_argument(\
        '--instance-type', metavar='STR', type=str, help='Amazon EC2 instance type to use for all nodes.')
    parser.add_argument(\
        '--instance-types', metavar='STR,STR,STR', type=str, help='Comma-separated list of EC2 instance types to use for MASTER, CORE and TASK nodes.')
    parser.add_argument(\
        '--instance-count', metavar='INT', type=int, help='Number of EC2 instances to use (1 MASTER, rest will be CORE).')
    parser.add_argument(\
        '--instance-counts', metavar='INT,INT,INT', type=str, help='Comma-separated list of number of EC2 instances to use for MASTER, CODE and TASK nodes.')
    parser.add_argument(\
        '--bid-price', metavar='FLOAT', type=float, help='EC2 spot instance bid price.  If this is specified, TASK nodes will be spot instances.')
