"""
rail-rna_config.py
Part of Rail-RNA

Contains classes that generate JSON configuration files for Rail-RNA from
application-specific command-line parameters. These configuration files are
parsable by Dooplicity's hadoop_simulator.py and emr_runner.py.

Class structure is designed so only those arguments relevant to scripts
are included. For example, if a script "rail-rna_preprocess_local.py"
were written, add_args(parser) for that class would add to the parser only
those arguments relevant to preprocessing input data in local mode.

The descriptions of command-line arguments contained here assume that a
calling script has the command-line options "--local", "--cloud",
"--preprocess", "--align", and "--all".

TO PUT IN rail-rna.py:

modes = set(['local', 'cloud'])
        # Implement Hadoop mode after paper submission
        if mode not in modes:
            self.errors.append('Mode ("--mode") must be one of '
                               '{{"local", "cloud"}}, but {0} was '
                               'entered.'.format(mode))
        self.mode = mode
        job_flows = set(['preprocess', 'align', 'all'])
        if job_flow not in job_flows:
            self.errors.append('Job flow ("--job-flows") must be one of '
                               '{{"preprocess", "align", "all"}}, but {0} was '
                               'entered.'.format(mode))
"""

import os

base_path = os.path.abspath(
                    os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
                )
utils_path = os.path.join(base_path, 'rna/utils')
site.addsitedir(utils_path)
site.addsitedir(base_path)
import dooplicity.ansibles as ab
from filemover import FileMover
import tempfile
import shutil

def add_args(parser):
    """ Adds relevant arguments to an object of class argparse.ArgumentParser.

        No return value.
    """
    parser.add_argument(
            '-b', '--branding', type=str, required=False, default=None,
            help='Text file with heading to write to stdout when running job. '
                 'This is where the name of a software package or ASCII art '
                 'can go.'
        )

class RailRnaJson(object):
    """ Base class for constructing JSON from input parameters.

        Checks only those parameters that are common to all modes/job flows.
    """
    def __init__(self, manifest, output_dir,
            intermediate_dir='./intermediate', force=False,
            s3cmd_exe='s3cmd', s3cmd_config=None
        ):
        # At the very least, the JSON payload will be composed of steps
        self.json_payload = { 'Steps': [] }
        '''Store all errors uncovered in a list, then output. This prevents the
        user from having to rerun Rail-RNA to find what else is wrong with
        the command-line parameters.'''
        self.errors = []
        self.temp = []
        self.manifest = manifest
        self.output_dir = output_dir
        self.intermediate_dir = intermediate_dir
        self.s3cmd_exe = s3cmd_exe
        self.s3cmd_config = s3cmd_config
        self.checked_s3 = False
        self.file_mover = FileMover(s3cmd_exe=self.s3cmd_exe,
                                    s3cred=self.s3cmd_config)
        self.initialized = True

    def check_s3(self, reason=None):
        """ Checks for s3cmd and AWS credentials; sets up S3Ansible.

            In this script, S3 checking is performed as soon as it is found
            that S3 is needed. If anything is awry, a RuntimeError is raised
            _immediately_ (the standard behavior is to raise a RuntimeError
            only after errors are accumulated). A reason specifying where
            S3 credentials were first needed can also be provided.

            reason: string specifying where S3 credentials were first
                needed.

            No return value.
        """
        if not ab.which(s3cmd_exe):
            self.errors.append(('The s3cmd executable ("--s3cmd-exe") "{0}" '
                                'was not found. Make sure that the '
                                'executable is in PATH or specify the '
                                'location of the executable with '
                                '"--s3cmd-exe".').format(self.s3cmd_exe))
        if self.s3cmd_config is not None:
            self.s3_ansible = ab.S3Ansible(s3cmd_exe=self.s3cmd_exe,
                                            config=self.s3cmd_config)
        elif os.name == 'nt' and os.getenv('USERPROFILE'):
            s3cmd_config_file = os.path.join(os.getenv('USERPROFILE'),
                    'Application Data', 's3cmd.ini'
                )
            self.s3_ansible = ab.S3Ansible(s3cmd_exe=self.s3cmd_exe,
                                        config=self.s3cmd_config)
        elif os.getenv('HOME'):
            s3cmd_config_file = os.path.join(os.getenv('HOME'), '.s3cfg')
            self.s3_ansible = ab.S3Ansible()
        else:
            self.errors.append('A valid s3cmd configuration file was not '
                               'found. Check that s3cmd is installed '
                               'properly, or specify the location of the '
                               'configuration file with "--s3cmd-config".')
        if reason:
            raise RuntimeError(('\n'.join(['%d) %s' % (i, error) for i, error
                                in enumerate(self.errors)]) + 
                                '\n\nNote that s3cmd is needed because {0}. '
                                'If all dependence on S3 in the pipeline is\n'
                                'removed, s3cmd need not be '
                                'installed.').format(reason))
        else:
            raise RuntimeError('\n'.join(['%d) %s' % (i, error) for i, error
                                in enumerate(self.errors)]) + 
                                '\n\nIf all dependence on S3 in the pipeline '
                                'is removed, s3cmd need not be installed.')
        self.checked_s3 = True

    def add_args(parser):
        parser.add_argument(
            '--s3cmd-exe', type=str, required=False,
            default='s3cmd',
            help='s3cmd executable. If s3cmd is in PATH, this is simply ' \
                 '"s3cmd". If it\'s not, the full path to s3cmd should be ' \
                 'specified. If S3 is never used in a given job flow, ' \
                 'this parameter is inconsequential.'
        )
        parser.add_argument(
            '--s3cmd-config', type=str, required=False,
            default=None,
            help='Path to s3cmd configuration file. If this is left ' \
                 'unspecified, the default configuration is used.'
        )

class RailRnaJsonLocal(RailRnaJson):
    """ Base class for constructing local-mode JSON from input parameters.

        Subsumes only those parameters relevant to local mode. Checks for files
        on local filesystem.
    """
    def __init__(self, manifest, output_dir, intermediate_dir='./intermediate',
                    force=False, s3cmd_exe='s3cmd', s3cmd_config=None,
                    num_processes=1, keep_intermediates=False):
        try:
            if self.initialized:
                pass
        except NameError:
            # Must initialize superclass
            RailRnaJson.__init__(manifest, output_dir,
                                    intermediate_dir=intermediate_dir,
                                    force=force, s3cmd_exe=s3cmd_exe,
                                    s3cmd_config=s3cmd_config)
            self.initialized = True
        output_dir_url = ab.Url(output_dir)
        if not (output_dir_url.is_local or output_dir_url.is_s3):
            self.errors.append(('Output directory must be local or on S3 '
                                'when running Rail-RNA in local ("--local") '
                                'mode, but {0} was entered.').format(
                                        output_dir
                                    ))
        if output_dir_url.is_s3:
            if not self.checked_s3:
                RailRnaJson.check_s3(reason='the output directory is on S3')
            # Create S3 bucket if it doesn't exist
            self.s3_ansible.create_bucket(ab.bucket_from_url(output_dir))
        if not ab.Url(intermediate_dir).is_local:
            self.errors.append(('Intermediate directory must be local '
                                'when running Rail-RNA in local ("--local") '
                                'mode, but {0} was entered.').format(
                                        intermediate_dir
                                    ))
        # Check manifest; download it if necessary
        manifest_url = ab.Url(self.manifest)
        if manifest_url.is_s3 and not self.checked_s3:
            RailRnaJson.check_s3(reason='the manifest file is on S3')
        if not file_mover.exists(manifest_url.to_url()):
            self.errors.append(('Manifest file ("--manifest") {0} '
                                'does not exist. Check the URL and '
                                'try again.').format(manifest))
        else:
            self.manifest_dir = tempfile.mkdtemp()
            if not manifest_url.is_local:
                self.manifest = os.path.join(self.manifest_dir, 'MANIFEST')
                file_mover.get(manifest_url, dest=self.manifest)
            else:
                self.manifest = manifest
            files_to_check = []
            with open(manifest) as manifest_stream:
                for line in manifest_stream:
                    tokens = line.strip().split('\t')
                    if len(tokens) == 5:
                        files_to_check.extend([tokens[0], tokens[2]])
                    elif len(tokens) == 3:
                        files_to_check.append(tokens[0])
                    else:
                        self.errors.append(('The following line from the '
                                            'manifest file {0} '
                                            'has an invalid number of '
                                            'tokens:\n{1}'
                                            ).format(
                                                    manifest,
                                                    line
                                                ))
            if files_to_check:
                for filename in files_to_check:
                    filename_url = ab.Url(filename)
                    if not file_mover.exists(filename_url):
                        self.errors.append(('The file {0} from the manifest '
                                            'file {1} does not exist. Check '
                                            'the URL and try again.').format(
                                                                filename,
                                                                manifest
                                                            ))
            else:
                self.errors.append(('Manifest file ("--manifest") {0} '
                                    'has no valid lines.').format(
                                                                manifest
                                                            ))
        from multiprocessing import cpu_count
        if num_processes:
            if not (isinstance(num_processes, int)
                                    and num_processes >= 1):
                self.errors.append('Number of processes ("--num-processes") '
                                   'must be an integer >= 1, '
                                   'but {0} was entered.'.format(
                                                    num_processes
                                                ))
            else:
                self.num_processes = num_processes
        else:
            try:
                self.num_processes = cpu_count()
            except NotImplementedError:
                self.num_processes = 1
            if self.num_processes != 1:
                self.num_processes -= 1
        self.keep_intermediates = keep_intermediates

    def add_args(parser):
        """ Adds parameter descriptions relevant to local mode to an object
            of class argparse.ArgumentParser.

            No return value.
        """
        parser.add_argument(
            '-m', '--manifest', type=str, required=True,
            help='Myrna-style manifest file listing sample FASTAs and/or ' \
                 'FASTQs to analyze. When running Rail-RNA in local ' \
                 '("--local") mode, this file must be on the local ' \
                 'filesystem. Each line of the file is ' \
                 'is in one of two formats: ' \
                 '\n\n(for unpaired input reads)\n' \
                 '<URL><TAB><MD5 checksum or "0" if not included><TAB>' \
                 '<sample label>' \
                 '\n\n(for paired-end input reads)\n' \
                 '<URL 1><TAB><MD5 checksum or "0" if not included><TAB>' \
                 '<URL 2><TAB><MD5 checksum or "0" if not included><TAB>' \
                 '<sample label>.'
        )
        parser.add_argument(
            '-o', '--output', type=str, required=False,
            default='./rail-rna_out',
            help='Output directory. This directory is not overwritten ' \
                 'unless the force overwrite ("--force") parameter is also ' \
                 'invoked.'
        )
        parser.add_argument(
            '--intermediate', type=str, required=False,
            default='./rail-rna_intermediate',
            help='Directory in which to store intermediate files, which ' \
                 'may be useful for debugging. Invoke ' \
                 '"--keep-intermediates" to prevent deletion of ' \
                 'intermediate files after a job flow is completed.'
        )
        parser.add_argument(
            '--num-processes', type=int, required=True,
            default=None,
            help='Number of processes to run simultaneously. This defaults to '
                 'the number of cores on the machine less 1 if more than one '
                 'core is available, or simply 1 if the program could not '
                 'determine the number of available cores.'
        )
        parser.add_argument(
            '--keep-intermediates', action='store_const', const=True,
            default=False,
            help='Keeps intermediate files after a job flow is completed.'
        )

class RailRnaJsonPreprocess(RailRnaJson):
    """ Base class for constructing preprocess-mode JSON from input parameters.
    """
    def __init__():
        super(RailRnaJsonPreprocess, self).__init__()

class RailRnaJsonAlign(RailRnaJson):

class RailRnaJsonAll(RailRnaJsonPreprocess, RailRnaJsonAlign):

class RailRnaJsonPreprocessLocal(RailRnaJson, RailRnaJsonPreprocess,
                                                RailRnaJsonLocal):

class RailRnaJsonAlignLocal(RailRnaJsonAlign, RailRnaJsonLocal):

class RailRnaJsonAllLocal(RailRnaJsonAll, RailRnaJsonLocal)

class RailRnaJsonPreprocessCloud(RailRnaJsonPreprocess, RailRnaJsonAws):


class RailRnaJson:
    """ Base class for constructing JSON from input parameters.

        Checks only those parameters that are common to all modes/job flows.
    """
    def __init__(self, manifest, output_dir, input_dir=None,
        intermediate_dir='intermediate', mode='local',
        job_flow='all', region='us-east-1', 
        log_uri=None, ami_version='2.4.2', visible_to_all_users=False, tags=[],
        name='Rail-RNA Job Flow',
        action_on_failure='TERMINATE_JOB_FLOW',
        hadoop_jar='/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar',
        num_processes=1, master_instance_count=1,
        master_instance_type='c1.xlarge',
        master_instance_bid_price=None, core_instance_count=1,
        core_instance_type=None, core_instance_bid_price=None,
        task_instance_count=0, task_instance_type=None,
        task_instance_bid_price=None, ec2_key_name=None, keep_alive=False,
        termination_protected=False, verbose=True,
        nucleotides_per_input=8000000, gzip_input=True, 
        bowtie_idx='genome', bowtie2_idx='genome', bowtie2_args='',
        genome_partition_length=5000, max_readlet_size=25, min_readlet_size=15,
        readlet_interval=4, cap_size_multiplier=1.2,
        max_intron_size=500000, min_intron_size=10, min_exon_size=9,
        motif_search_window_size=1000, motif_radius=5,
        normalize_percentile=0.75, output_bam_by_sample=False,
        bam_basename='alignments', bed_basename=''):

        # Initialize possible options
        self.instance_core_counts = {
            "m1.small"    : 1,
            "m1.large"    : 2,
            "m1.xlarge"   : 4,
            "c1.medium"   : 2,
            "c1.xlarge"   : 8,
            "m2.xlarge"   : 2,
            "m2.2xlarge"  : 4,
            "m2.4xlarge"  : 8,
            "cc1.4xlarge" : 8
        }

        self.instance_swap_allocations = {
            "m1.small"    : (2 *1024), #  1.7 GB
            "m1.large"    : (8 *1024), #  7.5 GB
            "m1.xlarge"   : (16*1024), # 15.0 GB
            "c1.medium"   : (2 *1024), #  1.7 GB
            "c1.xlarge"   : (8 *1024), #  7.0 GB
            "m2.xlarge"   : (16*1024), # 17.1 GB
            "m2.2xlarge"  : (16*1024), # 34.2 GB
            "m2.4xlarge"  : (16*1024), # 68.4 GB
            "cc1.4xlarge" : (16*1024)  # 23.0 GB
        }

        '''Not currently in use, but may become important if there are
        32- vs. 64-bit issues: self.instance_bits = {
            "m1.small"    : 32,
            "m1.large"    : 64,
            "m1.xlarge"   : 64,
            "c1.medium"   : 32,
            "c1.xlarge"   : 64,
            "m2.xlarge"   : 64,
            "m2.2xlarge"  : 64,
            "m2.4xlarge"  : 64,
            "cc1.4xlarge" : 64
        }'''

        # Perform all command-line parameter checking here
        modes = set(['local', 'cloud'])
        # Implement Hadoop mode after paper submission
        if mode not in modes:
            self.errors.append('Mode ("--mode") must be one of '
                               '{"local", "cloud"}, but {0} was '
                               'entered.'.format(mode))
        self.mode = mode
        job_flows = set(['preprocess', 'align', 'all'])
        if job_flow not in job_flows:
            self.errors.append('Job flow ("--job-flows") must be one of '
                               '{"preprocess", "align", "all"}, but {0} was '
                               'entered.'.format(mode))
        self.job_flow = job_flow
        actions_on_failure \
            = set(['TERMINATE_JOB_FLOW', 'CANCEL_AND_WAIT', 'CONTINUE',
                    'TERMINATE_CLUSTER'])
        self.manifest = manifest
        self.intermediate_dir = intermediate_dir
        self.output_dir = output_dir
        if action_on_failure not in actions_on_failure:
            self.errors.append('Action on failure ("--action-on-failure") '
                               'must be one of {"TERMINATE_JOB_FLOW", '
                               '"CANCEL_AND_WAIT", "CONTINUE", '
                               '"TERMINATE_CLUSTER"}, but '
                               '{0} was entered.'.format(
                                                action_on_failure
                                            ))
        self.action_on_failure = action_on_failure
        self.hadoop_jar = hadoop_jar
        if not (isinstance(num_processes, int)
                and num_processes >= 1):
            self.errors.append('Number of processes ("--num-processes") must '
                               'be an integer >= 1, '
                               'but {0} was entered.'.format(
                                                num_processes
                                            ))
        self.tasks_per_reducer = tasks_per_reducer
        self.reducer_count = reducer_count
        instance_type_message = ('Instance type ("--instance-type") must be '
                                 'in the set {"m1.small", "m1.large", '
                                 '"m1.xlarge", "c1.medium", "c1.xlarge", '
                                 '"m2.xlarge", "m2.2xlarge", "m2.4xlarge", '
                                 '"cc1.4xlarge"}, but {0} was entered.')
        if master_instance_type not in instance_core_counts:
            self.errors.append(('Master instance type '
                               '("--master-instance-type") not valid. %s')
                                % instance_type_message.format(
                                                        master_instance_type
                                                    ))
        self.master_instance_type = master_instance_type
        if core_instance_type is None:
            self.core_instance_type = self.master_instance_type
        else:
            if core_instance_type not in instance_core_counts:
                self.errors.append(('Core instance type '
                                    '("--core-instance-type") not valid. %s')
                                    % instance_type_message.format(
                                                        core_instance_type
                                                    ))
            self.core_instance_type = core_instance_type
        if task_instance_type is None:
            self.task_instance_type = self.master_instance_type
        else:
            if task_instance_type not in instance_core_counts:
                self.errors.append(('Task instance type '
                                    '("--task-instance-type") not valid. %s')
                                    % instance_type_message.format(
                                                        task_instance_type
                                                    ))
            self.task_instance_type = task_instance_type
        if master_instance_bid_price is None:
            self.spot_master = False
        else:
            if not (isinstance(master_instance bid_price, float) 
                    and master_instance_bid_price > 0):
                self.errors.append('Spot instance bid price for master nodes '
                                   '(--master-instance-bid-price) must be '
                                   '> 0, but {0} was entered.'.format(
                                                    master_instance_bid_price
                                                ))
            self.spot_master = True
            self.master_instance_bid_price = master_instance_bid_price
        if core_instance_bid_price is None:
            self.spot_core = False
        else:
            if not (isinstance(core_instance bid_price, float) 
                    and core_instance_bid_price > 0):
                self.errors.append('Spot instance bid price for core nodes '
                                   '(--core-instance-bid-price) must be '
                                   '> 0, but {0} was entered.'.format(
                                                    core_instance_bid_price
                                                ))
            self.spot_core = True
            self.core_instance_bid_price = core_instance_bid_price
        if task_instance_bid_price is None:
            self.spot_task = False
        else:
            if not (isinstance(task_instance bid_price, float) 
                    and task_instance_bid_price > 0):
                self.errors.append('Spot instance bid price for task nodes '
                                   '(--task-instance-bid-price) must be '
                                   '> 0, but {0} was entered.'.format(
                                                    task_instance_bid_price
                                                ))
            self.spot_task = True
            self.task_instance_bid_price = task_instance_bid_price
        if not (isinstance(master_instance_count, int)
                and master_instance_count >= 1):
            self.errors.append('Master instance count '
                               '("--master-instance-count") must be an '
                               'integer >= 1, but {0} was entered.'.format(
                                                    master_instance_count
                                                ))
        self.master_instance_count = master_instance_count
        if not (isinstance(core_instance_count, int)
                and core_instance_count >= 0):
            self.errors.append('Core instance count '
                               '("--core-instance-count") must be an '
                               'integer >= 1, but {0} was entered.'.format(
                                                    core_instance_count
                                                ))
        self.core_instance_count = core_instance_count
        if not (isinstance(task_instance_count, int)
                and task_instance_count >= 0):
            self.errors.append('Task instance count '
                               '("--task-instance-count") must be an '
                               'integer >= 1, but {0} was entered.'.format(
                                                    task_instance_count
                                                ))
        self.task_instance_count = task_instance_count
        if self.core_instance_count > 0:
            self.swap_allocation \
                = self.instance_swap_allocations[self.core_instance_type]
        else:
            self.swap_allocation \
                = self.instance_swap_allocations[self.master_instance_type]
        self.ec2_key_name = ec2_key_name
        self.keep_alive = keep_alive
        self.termination_protection = termination_protection
        self.verbose = verbose
        if not (isinstance(nucleotides_per_input, int) and
                nucleotides_per_input > 0):
            self.errors.append('Nucleotides per input '
                               '(--nucleotides-per-input) must be an integer '
                               '> 0, but {0} was entered.'.format(
                                                        nucleotides_per_input
                                                       ))
        self.nucleotides_per_input = nucleotides_per_input
        if not (isinstance(genome_partition_length, int) and
                genome_partition_length > 0):
            self.errors.append('Genome partition length '
                               '(--genome-partition-length) must be an '
                               'integer > 0, but {0} was entered.'.format(
                                                        genome_partition_length
                                                    ))
        self.genome_partition_length = genome_partition_length
        if not (isinstance(min_readlet_size, int) and min_readlet_size > 0):
            self.errors.append('Minimum readlet size (--min-readlet-size) '
                               'must be an integer > 0, but '
                               '{0} was entered.'.format(min_readlet_size))
        self.min_readlet_size = min_readlet_size
        if not (isinstance(max_readlet_size, int) and max_readlet_size
                >= min_readlet_size):
            self.errors.append('Maximum readlet size (--max-readlet-size) '
                               'must be an integer >= minimum readlet size '
                               '(--min-readlet-size) = '
                               '{0}, but {1} was entered.'.format(
                                                    self.min_readlet_size,
                                                    max_readlet_size
                                                ))
        self.max_readlet_size = max_readlet_size
        if not (isinstance(readlet_interval, int) and readlet_interval
                > 0):
            self.errors.append('Readlet interval (--readlet-interval) '
                               'must be an integer > 0, '
                               'but {0} was entered.'.format(
                                                    readlet_interval
                                                ))
        self.readlet_interval = readlet_interval
        if not (isinstance(cap_size_multiplier, float) and cap_size_multiplier
                > 1):
            self.errors.append('Cap size multiplier (--cap-size-multiplier) '
                               'must be > 1, '
                               'but {0} was entered.'.format(
                                                    cap_size_multiplier
                                                ))
        self.cap_size_multiplier = cap_size_multiplier
        if not (isinstance(min_intron_size, int) and min_intron_size > 0):
            self.errors.append('Minimum intron size (--min-intron-size) '
                               'must be an integer > 0, but '
                               '{0} was entered.'.format(min_intron_size))
        self.min_intron_size = min_intron_size
        if not (isinstance(max_intron_size, int) and max_intron_size
                >= min_intron_size):
            self.errors.append('Maximum intron size (--max-intron-size) '
                               'must be an integer >= minimum intron size '
                               '(--min-readlet-size) = '
                               '{0}, but {1} was entered.'.format(
                                                    self.min_intron_size,
                                                    max_intron_size
                                                ))
        self.max_intron_size = max_intron_size
        if not (isinstance(min_exon_size, int) and min_exon_size > 0):
            self.errors.append('Minimum exon size (--min-exon-size) '
                               'must be an integer > 0, but '
                               '{0} was entered.'.format(min_exon_size))
        self.min_exon_size = min_exon_size
        if not (isinstance(motif_search_window_size, int) and 
                    motif_search_window_size >= 0):
            self.errors.append('Motif search window size '
                               '(--motif-search-window-size) must be an '
                               'integer > 0, but {0} was entered.'.format(
                                                    motif_search_window_size
                                                ))
        if not (isinstance(motif_radius, int) and
                    motif_radius >= 0):
            self.errors.append('Motif radius (--motif-radius) must be an '
                               'integer >= 0, but {0} was entered.'.format(
                                                    motif_radius
                                                ))
        if not (isinstance(normalize_percentile, float) and
                    0 <= normalize_percentile <= 1):
            self.errors.append('Normalization percentile '
                               '(--normalize-percentile) must on the '
                               'interval [0, 1], but {0} was entered'.format(
                                                    normalize_percentile
                                                ))
        self.bam_basename = bam_basename
        self.bed_basename = bed_basename
        self.preprocess_bootstrap_json = \
 """[
        {{
            "Name": "Install PyPy",
            "ScriptBootstrapAction": {{
                "Args": [
                    "s3://rail-emr/bin/pypy-2.2.1-linux_x86_64-portable.tar.bz2"
                ],
                "Path": "s3://rail-emr/bootstrap/install-pypy.sh"
            }}
        }},
        {{
            "Name": "Install Rail-RNA",
            "ScriptBootstrapAction": {{
                "Args": [
                    "s3://rail-emr/bin/rail-rna-0.1.0.tar.gz",
                    "/mnt"
                ],
                "Path": "s3://rail-emr/bootstrap/install-rail.sh"
            }}
        }},
        {{
            "Name": "Install manifest file",
            "ScriptBootstrapAction": {{
                "Args": [
                    "s3://rail-experiments/geuvadis_abbreviated/manifest.20samples",
                    "/mnt",
                    "MANIFEST"
                ],
                "Path": "s3://rail-emr/bootstrap/s3cmd_s3.sh"
            }}
        }},
        {{
            "Name": "Add swap space",
            "ScriptBootstrapAction": {{
                "Args": [
                    "{swap_allocation}"
                ],
                "Path": "s3://elasticmapreduce/bootstrap-actions/add-swap"
            }}
        }},
        {{
            "Name": "Configure Hadoop",
            "ScriptBootstrapAction": {{
                "Args": [
                    "-s", "mapred.job.reuse.jvm.num.tasks=1",
                    "-s", "mapred.tasktracker.reduce.tasks.maximum={instance_core_count}",
                    "-s", "mapred.tasktracker.map.tasks.maximum={instance_core_count}",
                    "-m", "mapred.map.tasks.speculative.execution=false",
                    "-m", "mapred.reduce.tasks.speculative.execution=false"
                ],
                "Path": "s3://elasticmapreduce/bootstrap-actions/configure-hadoop"
            }}
        }}
    ]""".format(swap_allocation=self.swap_allocation,
                instance_core_count=(
                      self.instance_core_counts[self.core_instance_type]
                      if core_instance_count > 0 
                      else self.instance_core_counts[self.master_instance_type]
                    )
                )
        self.process_bootstrap_json = \
 """[
        {{
            "Name": "Install PyPy",
            "ScriptBootstrapAction": {{
                "Args": [
                    "s3://rail-emr/bin/pypy-2.2.1-linux_x86_64-portable.tar.bz2"
                ],
                "Path": "s3://rail-emr/bootstrap/install-pypy.sh"
            }}
        }},
        {{
            "Name": "Install Bowtie 1",
            "ScriptBootstrapAction": {{
                "Args": [],
                "Path": "s3://rail-emr/bootstrap/install-bowtie.sh"
            }}
        }},
        {{
            "Name": "Install Bowtie 2",
            "ScriptBootstrapAction": {{
                "Args": [],
                "Path": "s3://rail-emr/bootstrap/install-bowtie2.sh"
            }}
        }},
        {{
            "Name": "Install bedToBigBed",
            "ScriptBootstrapAction": {{
                "Args": [
                    "/mnt/bin"
                ],
                "Path": "s3://rail-emr/bootstrap/install-kenttools.sh"
            }}
        }},
        {{
            "Name": "Install SAMtools",
            "ScriptBootstrapAction": {{
                "Args": [],
                "Path": "s3://rail-emr/bootstrap/install-samtools.sh"
            }}
        }},
        {{
            "Name": "Install Rail-RNA",
            "ScriptBootstrapAction": {{
                "Args": [
                    "s3://rail-emr/bin/rail-rna-0.1.0.tar.gz",
                    "/mnt"
                ],
                "Path": "s3://rail-emr/bootstrap/install-rail.sh"
            }}
        }},
        {{
            "Name": "Install Bowtie indexes",
            "ScriptBootstrapAction": {{
                "Args": [
                    "/mnt",
                    "s3://rail-emr/index/hg19_UCSC.tar.gz"
                ],
                "Path": "s3://rail-emr/bootstrap/s3cmd_s3_tarball.sh"
            }}
        }},
        {{
            "Name": "Install manifest file",
            "ScriptBootstrapAction": {{
                "Args": [
                    "s3://rail-experiments/geuvadis_abbreviated/manifest.20samples",
                    "/mnt",
                    "MANIFEST"
                ],
                "Path": "s3://rail-emr/bootstrap/s3cmd_s3.sh"
            }}
        }},
        {{
            "Name": "Add swap space",
            "ScriptBootstrapAction": {{
                "Args": [
                    "{swap_allocation}"
                ],
                "Path": "s3://elasticmapreduce/bootstrap-actions/add-swap"
            }}
        }},
        {{
            "Name": "Configure Hadoop",
            "ScriptBootstrapAction": {{
                "Args": [
                    "-s", "mapred.job.reuse.jvm.num.tasks=1",
                    "-s", "mapred.tasktracker.reduce.tasks.maximum={instance_core_count}",
                    "-s", "mapred.tasktracker.map.tasks.maximum={instance_core_count}",
                    "-m", "mapred.map.tasks.speculative.execution=false",
                    "-m", "mapred.reduce.tasks.speculative.execution=false"
                ],
                "Path": "s3://elasticmapreduce/bootstrap-actions/configure-hadoop"
            }}
        }}
    ]""".format(swap_allocation=self.swap_allocation,
                instance_core_count=(
                      self.instance_core_counts[self.core_instance_type]
                      if core_instance_count > 0 
                      else self.instance_core_counts[self.master_instance_type]
                    )
                )
        self.instance_json = \
 """{
        {ec2_key_name}
        "HadoopVersion": "1.0.3",
        "InstanceGroups": [
            {
                {master_spot}
                "InstanceCount": {master_instance_count},
                "InstanceRole": "MASTER",
                "InstanceType": "{master_instance_type}",
                "Name": "Master Instance Group"
            },
            {
                {core_spot}
                "InstanceCount": {core_instance_count},
                "InstanceRole": "CORE",
                "InstanceType": "{core_instance_type}",
                "Name": "Core Instance Group"
            },
            {
                {task_spot}
                "InstanceCount": {task_instance_count},
                "InstanceRole": "TASK",
                "InstanceType": "{task_instance_type}",
                "Name": "Task Instance Group"
            }
        ],
        "KeepJobFlowAliveWhenNoSteps": "false",
        "TerminationProtected": "{termination_protected}"
    }""".format(
                ec2_key_name=(
                        ('"Ec2KeyName": "%s",' % self.ec2_key_name)
                         if self.ec2_key_name is not None else ''
                    ),
                master_spot=(
                        (('"BidPrice": "{%0.4f}",' 
                            % self.master_instance_bid_price)
                          '"Market": "SPOT",')
                        if self.master_instance_bid_price is not None else
                        '"Market" : "ON_DEMAND",'
                    ),
                core_spot=(
                        (('"BidPrice": "{%0.4f}",' 
                            % self.core_instance_bid_price)
                          '"Market": "SPOT",')
                        if self.core_instance_bid_price is not None else
                        '"Market" : "ON_DEMAND",'
                    ),
                task_spot=(
                        (('"BidPrice": "{%0.4f}",' 
                            % self.task_instance_bid_price)
                          '"Market": "SPOT",')
                        if self.task_instance_bid_price is not None else
                        '"Market" : "ON_DEMAND",'
                    )
                master_instance_count=self.master_instance_count,
                core_instance_count=self.core_instance_count,
                task_instance_count=self.task_instance_count,
                master_instance_type=self.master_instance_type,
                core_instance_type=self.core_instance_type,
                task_instance_type=self.task_instance_type,
                termination_protected=('true' if self.termination_protected
                                        else 'false')
            ),
        self.hadoop_debugging_step = \
    """{
            "Name": "Set up Hadoop Debugging"
            "ActionOnFailure": {action_on_failure},
            "HadoopJarStep": {
                "Args": [
                    "s3://us-east-1.elasticmapreduce/libs/state-pusher/0.1/fetch"
                ],
                "Jar": "s3://us-east-1.elasticmapreduce/libs/script-runner/script-runner.jar"
            },
        }""".format(action_on_failure=self.action_on_failure)
        self.cloud_preprocess_step = \
    """{
            "Name": "Preprocess input reads and store them on S3"
            "ActionOnFailure": {action_on_failure},
            "HadoopJarStep": {
                "Args": [
                    "-D", "mapred.reduce.tasks=0",
                    "-input", "{manifest}",
                    "-output", "{preprocess_dir}",
                    "-mapper", "pypy /mnt/src/rail-rna/preprocess.py --nucs-per-file={nucleotides_per_input} {gzip_output} --push={upload_dir} --ignore-first-token",
                    "-reducer", "cat",
                    "-inputformat", "org.apache.hadoop.mapred.lib.NLineInputFormat"
                ],
                "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar"
            },
        }""".format(action_on_failure=self.action_on_failure,
                    manifest=self.manifest
                    preprocess_dir=os.path.join(self.intermediate_dir, 
                                                'preprocess')
                    nucleotides_per_input=self.nucleotides_per_input,
                    gzip_output=('--gzip-output' if self.gzip_output else ''),
                    upload_dir=(
                        os.path.join(self.intermediate_dir, 'preprocess/push')
                        if self.job_flow != 'preprocess'
                        else os.path.join(self.output_dir, '')
                    )
                )
        self.cloud_align_steps = \
 """[
        {
            "Name" : "Align reads to genome"
            "ActionOnFailure": {action_on_failure},
            "HadoopJarStep": {
                "Args": [
                    "-libjars", "/mnt/lib/multiplefiles.jar",
                    "-D", "mapred.reduce.tasks=0",
                    "-input", "{align_input}",
                    "-output", "{align_output}",
                    "-mapper", "pypy /mnt/src/rail-rna/align_reads.py --bowtie-idx=/mnt/index/genome --bowtie2-idx=/mnt/index/genome --bowtie2-exe=bowtie2 --exon-differentials --partition-length {genome_partition_length} --manifest=/mnt/MANIFEST {verbose} -- --local",
                    "-outputformat", "edu.jhu.cs.MultipleOutputFormat"
                ],
                "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar"
            }
        },
        {
            "Name": "Aggregate duplicate read sequences to reduce realignment burden",
            "ActionOnFailure": {action_on_failure},
            "HadoopJarStep": {
                "Args": [
                    "-D", "mapred.reduce.tasks={combine_sequences_task_count}",
                    "-D", "mapred.text.key.partitioner.options=-k1,1",
                    "-D", "stream.num.map.output.key.fields=1",
                    "-partitioner", "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
                    "-input", "s3://rail-experiments/geuvadis_again_intermediate/align_reads/readletize",
                    "-output","s3://rail-experiments/geuvadis_again_intermediate/combine_sequences",
                    "-mapper",
                    "cat",
                    "-reducer",
                    "pypy /mnt/src/rail-rna/sum.py --type 3 --value-count 2"
                ],
                "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar"
            }
        },
        {
            "ActionOnFailure": {action_on_failure},
            "HadoopJarStep": {
                "Args": [
                    "-D",
                    "mapred.reduce.tasks=1600",
                    "-D",
                    "mapred.text.key.partitioner.options=-k1,1",
                    "-D",
                    "stream.num.map.output.key.fields=1",
                    "-partitioner",
                    "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
                    "-input",
                    "s3://rail-experiments/geuvadis_again_intermediate/combine_sequences/",
                    "-output",
                    "s3://rail-experiments/geuvadis_again_intermediate/readletize",
                    "-mapper",
                    "cat",
                    "-reducer",
                    "pypy /mnt/src/rail-rna/readletize.py --max-readlet-size 23 --readlet-interval 4 --capping-multiplier 1.200000"
                ],
                "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar"
            },
            "Name": "Readletize"
        },
        {
            "ActionOnFailure": {action_on_failure},
            "HadoopJarStep": {
                "Args": [
                    "-D",
                    "mapred.reduce.tasks=1600",
                    "-D",
                    "mapred.text.key.partitioner.options=-k1,1",
                    "-D",
                    "stream.num.map.output.key.fields=1",
                    "-partitioner",
                    "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
                    "-input",
                    "s3://rail-experiments/geuvadis_again_intermediate/readletize/",
                    "-output",
                    "s3://rail-experiments/geuvadis_again_intermediate/combine_subsequences",
                    "-mapper",
                    "cat",
                    "-reducer",
                    "pypy /mnt/src/rail-rna/sum.py --type 3"
                ],
                "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar"
            },
            "Name": "CombineSubsequences"
        },
        {
            "ActionOnFailure": {action_on_failure},
            "HadoopJarStep": {
                "Args": [
                    "-D",
                    "mapred.reduce.tasks=1600",
                    "-D",
                    "mapred.text.key.partitioner.options=-k1,1",
                    "-D",
                    "stream.num.map.output.key.fields=1",
                    "-partitioner",
                    "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
                    "-input",
                    "s3://rail-experiments/geuvadis_again_intermediate/combine_subsequences/",
                    "-output",
                    "s3://rail-experiments/geuvadis_again_intermediate/align_readlets",
                    "-mapper",
                    "cat",
                    "-reducer",
                    "pypy /mnt/src/rail-rna/align_readlets.py --bowtie-idx=/mnt/index/genome --bowtie-exe=bowtie --verbose -- -t --sam-nohead --startverbose -v 0 -a -m 80"
                ],
                "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar"
            },
            "Name": "AlignReadlets"
        },
        {
            "ActionOnFailure": {action_on_failure},
            "HadoopJarStep": {
                "Args": [
                    "-D",
                    "mapred.reduce.tasks=1600",
                    "-D",
                    "mapred.text.key.partitioner.options=-k1,1",
                    "-D",
                    "stream.num.map.output.key.fields=1",
                    "-partitioner",
                    "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
                    "-input",
                    "s3://rail-experiments/geuvadis_again_intermediate/align_readlets/",
                    "-output",
                    "s3://rail-experiments/geuvadis_again_intermediate/intron_search",
                    "-mapper",
                    "cat",
                    "-reducer",
                    "pypy /mnt/src/rail-rna/intron_search.py --bowtie-idx=/mnt/index/genome --partition-length 5000 --max-intron-size 500000 --min-intron-size 10 --min-exon-size 9 --search-window-size 1000 --motif-radius 5 --verbose"
                ],
                "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar"
            },
            "Name": "IntronSearch"
        },
        {
            "ActionOnFailure": {action_on_failure},
            "HadoopJarStep": {
                "Args": [
                    "-D",
                    "mapred.reduce.tasks=400",
                    "-D",
                    "mapred.text.key.partitioner.options=-k1,2",
                    "-D",
                    "stream.num.map.output.key.fields=4",
                    "-partitioner",
                    "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
                    "-input",
                    "s3://rail-experiments/geuvadis_again_intermediate/intron_search/",
                    "-output",
                    "s3://rail-experiments/geuvadis_again_intermediate/intron_config",
                    "-mapper",
                    "cat",
                    "-reducer",
                    "pypy /mnt/src/rail-rna/intron_config.py --readlet-size 23 --verbose"
                ],
                "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar"
            },
            "Name": "IntronConfig"
        },
        {
            "ActionOnFailure": {action_on_failure},
            "HadoopJarStep": {
                "Args": [
                    "-D",
                    "mapred.reduce.tasks=3200",
                    "-D",
                    "mapred.text.key.partitioner.options=-k1,4",
                    "-D",
                    "stream.num.map.output.key.fields=4",
                    "-partitioner",
                    "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
                    "-input",
                    "s3://rail-experiments/geuvadis_again_intermediate/intron_config/",
                    "-output",
                    "s3://rail-experiments/geuvadis_again_intermediate/intron_fasta",
                    "-mapper",
                    "cat",
                    "-reducer",
                    "pypy /mnt/src/rail-rna/intron_fasta.py --verbose --bowtie-idx=/mnt/index/genome"
                ],
                "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar"
            },
            "Name": "IntronFasta"
        },
        {
            "ActionOnFailure": {action_on_failure},
            "HadoopJarStep": {
                "Args": [
                    "-D",
                    "mapred.reduce.tasks=1",
                    "-D",
                    "mapred.text.key.partitioner.options=-k1,1",
                    "-D",
                    "stream.num.map.output.key.fields=1",
                    "-partitioner",
                    "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
                    "-input",
                    "s3://rail-experiments/geuvadis_again_intermediate/intron_fasta/",
                    "-output",
                    "s3://rail-experiments/geuvadis_again_intermediate/intron_index",
                    "-mapper",
                    "cat",
                    "-reducer",
                    "pypy /mnt/src/rail-rna/intron_index.py --bowtie-build-exe=bowtie-build --bowtie-idx=/mnt/index/genome --out=S3://rail-experiments/geuvadis_again/index --keep-alive"
                ],
                "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar"
            },
            "Name": "IntronIndex"
        },
        {
            "ActionOnFailure": {action_on_failure},
            "HadoopJarStep": {
                "Args": [
                    "-D",
                    "mapred.reduce.tasks=1600",
                    "-D",
                    "mapred.text.key.partitioner.options=-k1,1",
                    "-D",
                    "stream.num.map.output.key.fields=1",
                    "-archives",
                    "s3n://rail-experiments/geuvadis_again/index/intron.tar.gz#intron",
                    "-partitioner",
                    "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
                    "-input",
                    "s3://rail-experiments/geuvadis_again_intermediate/combine_subsequences/",
                    "-output",
                    "s3://rail-experiments/geuvadis_again_intermediate/realign_readlets",
                    "-mapper",
                    "cat",
                    "-reducer",
                    "pypy /mnt/src/rail-rna/align_readlets.py --bowtie-idx=intron/intron --bowtie-exe=bowtie --verbose -- -t --sam-nohead --startverbose -v 0 -a -m 80"
                ],
                "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar"
            },
            "Name": "RealignReadlets"
        },
        {
            "ActionOnFailure": {action_on_failure},
            "HadoopJarStep": {
                "Args": [
                    "-D",
                    "mapred.reduce.tasks=1600",
                    "-D",
                    "mapred.text.key.partitioner.options=-k1,1",
                    "-D",
                    "stream.num.map.output.key.fields=1",
                    "-partitioner",
                    "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
                    "-input",
                    "s3://rail-experiments/geuvadis_again_intermediate/realign_readlets/,s3://rail-experiments/geuvadis_again_intermediate/align_readlets/",
                    "-output",
                    "s3://rail-experiments/geuvadis_again_intermediate/cointron_search",
                    "-mapper",
                    "cat",
                    "-reducer",
                    "pypy /mnt/src/rail-rna/cointron_search.py --verbose"
                ],
                "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar"
            },
            "Name": "CointronSearch"
        },
        {
            "ActionOnFailure": {action_on_failure},
            "HadoopJarStep": {
                "Args": [
                    "-D",
                    "mapred.reduce.tasks=3200",
                    "-D",
                    "mapred.text.key.partitioner.options=-k1,2",
                    "-D",
                    "stream.num.map.output.key.fields=2",
                    "-partitioner",
                    "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
                    "-input",
                    "s3://rail-experiments/geuvadis_again_intermediate/cointron_search/",
                    "-output",
                    "s3://rail-experiments/geuvadis_again_intermediate/cointron_fasta",
                    "-mapper",
                    "cat",
                    "-reducer",
                    "pypy /mnt/src/rail-rna/cointron_fasta.py --verbose --bowtie-idx=/mnt/index/genome"
                ],
                "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar"
            },
            "Name": "CointronFasta"
        },
        {
            "ActionOnFailure": {action_on_failure},
            "HadoopJarStep": {
                "Args": [
                    "-D",
                    "mapred.reduce.tasks=1600",
                    "-D",
                    "mapred.text.key.partitioner.options=-k1,1",
                    "-D",
                    "stream.num.map.output.key.fields=1",
                    "-libjars",
                    "/mnt/lib/multiplefiles.jar",
                    "-partitioner",
                    "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
                    "-input",
                    "s3://rail-experiments/geuvadis_again_intermediate/align_reads/unmapped,s3://rail-experiments/geuvadis_again_intermediate/cointron_fasta/",
                    "-output",
                    "s3://rail-experiments/geuvadis_again_intermediate/realign_reads",
                    "-mapper",
                    "cat",
                    "-reducer",
                    "pypy /mnt/src/rail-rna/realign_reads.py --original-idx=/mnt/index/genome --bowtie2-exe=bowtie2 --partition-length 5000 --exon-differentials --manifest=/mnt/MANIFEST --verbose --keep-alive -- --end-to-end",
                    "-outputformat",
                    "edu.jhu.cs.MultipleOutputFormat"
                ],
                "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar"
            },
            "Name": "RealignReads"
        },
        {
            "ActionOnFailure": {action_on_failure},
            "HadoopJarStep": {
                "Args": [
                    "-D",
                    "mapred.reduce.tasks=3200",
                    "-D",
                    "mapred.text.key.partitioner.options=-k1,3",
                    "-D",
                    "stream.num.map.output.key.fields=3",
                    "-partitioner",
                    "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
                    "-input",
                    "s3://rail-experiments/geuvadis_again_intermediate/align_reads/exon_diff,s3://rail-experiments/geuvadis_again_intermediate/realign_reads/exon_diff",
                    "-output",
                    "s3://rail-experiments/geuvadis_again_intermediate/collapse",
                    "-mapper",
                    "cat",
                    "-reducer",
                    "pypy /mnt/src/rail-rna/sum.py"
                ],
                "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar"
            },
            "Name": "Collapse"
        },
        {
            "ActionOnFailure": {action_on_failure},
            "HadoopJarStep": {
                "Args": [
                    "-D",
                    "mapred.reduce.tasks=3200",
                    "-D",
                    "mapred.text.key.partitioner.options=-k1,2",
                    "-D",
                    "stream.num.map.output.key.fields=3",
                    "-libjars",
                    "/mnt/lib/multiplefiles.jar",
                    "-partitioner",
                    "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
                    "-input",
                    "s3://rail-experiments/geuvadis_again_intermediate/collapse/",
                    "-output",
                    "s3://rail-experiments/geuvadis_again_intermediate/coverage_pre",
                    "-mapper",
                    "cat",
                    "-reducer",
                    "pypy /mnt/src/rail-rna/coverage_pre.py --bowtie-idx=/mnt/index/genome --partition-stats",
                    "-outputformat",
                    "edu.jhu.cs.MultipleOutputFormat"
                ],
                "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar"
            },
            "Name": "CoveragePre"
        },
        {
            "ActionOnFailure": {action_on_failure},
            "HadoopJarStep": {
                "Args": [
                    "-D",
                    "mapred.reduce.tasks=400",
                    "-D",
                    "mapred.text.key.partitioner.options=-k1,1",
                    "-D",
                    "stream.num.map.output.key.fields=3",
                    "-partitioner",
                    "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
                    "-input",
                    "s3://rail-experiments/geuvadis_again_intermediate/coverage_pre/coverage",
                    "-output",
                    "s3://rail-experiments/geuvadis_again_intermediate/coverage",
                    "-mapper",
                    "cat",
                    "-reducer",
                    "pypy /mnt/src/rail-rna/coverage.py --bowtie-idx=/mnt/index/genome --percentile 0.750000 --out=S3://rail-experiments/geuvadis_again/coverage --bigbed-exe=/mnt/bin/bedToBigBed --manifest=/mnt/MANIFEST --keep-alive --verbose"
                ],
                "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar"
            },
            "Name": "Coverage"
        },
        {
            "ActionOnFailure": {action_on_failure},
            "HadoopJarStep": {
                "Args": [
                    "-D",
                    "mapred.reduce.tasks=1",
                    "-D",
                    "mapred.text.key.partitioner.options=-k1,1",
                    "-D",
                    "stream.num.map.output.key.fields=2",
                    "-partitioner",
                    "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
                    "-input",
                    "s3://rail-experiments/geuvadis_again_intermediate/coverage",
                    "-output",
                    "s3://rail-experiments/geuvadis_again_intermediate/coverage_post",
                    "-mapper",
                    "cat",
                    "-reducer",
                    "pypy /mnt/src/rail-rna/coverage_post.py --out=S3://rail-experiments/geuvadis_again/normalize --manifest=/mnt/MANIFEST"
                ],
                "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar"
            },
            "Name": "CoveragePost"
        },
        {
            "ActionOnFailure": {action_on_failure},
            "HadoopJarStep": {
                "Args": [
                    "-D",
                    "mapred.reduce.tasks=3200",
                    "-D",
                    "mapred.text.key.partitioner.options=-k1,6",
                    "-D",
                    "stream.num.map.output.key.fields=6",
                    "-partitioner",
                    "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
                    "-input",
                    "s3://rail-experiments/geuvadis_again_intermediate/realign_reads/bed,s3://rail-experiments/geuvadis_again_intermediate/align_reads/bed",
                    "-output",
                    "s3://rail-experiments/geuvadis_again_intermediate/bed_pre",
                    "-mapper",
                    "cat",
                    "-reducer",
                    "pypy /mnt/src/rail-rna/bed_pre.py"
                ],
                "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar"
            },
            "Name": "BedPre"
        },
        {
            "ActionOnFailure": {action_on_failure},
            "HadoopJarStep": {
                "Args": [
                    "-D",
                    "mapred.reduce.tasks=400",
                    "-D",
                    "mapred.text.key.partitioner.options=-k1,2",
                    "-D",
                    "stream.num.map.output.key.fields=4",
                    "-partitioner",
                    "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
                    "-input",
                    "s3://rail-experiments/geuvadis_again_intermediate/bed_pre/",
                    "-output",
                    "s3://rail-experiments/geuvadis_again_intermediate/bed",
                    "-mapper",
                    "cat",
                    "-reducer",
                    "pypy /mnt/src/rail-rna/bed.py --bowtie-idx=/mnt/index/genome --out=S3://rail-experiments/geuvadis_again/bed --manifest=/mnt/MANIFEST --bed-basename="
                ],
                "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar"
            },
            "Name": "Bed"
        },
        {
            "ActionOnFailure": {action_on_failure},
            "HadoopJarStep": {
                "Args": [
                    "-D",
                    "mapred.reduce.tasks=400",
                    "-D",
                    "mapred.text.key.partitioner.options=-k1,1",
                    "-D",
                    "stream.num.map.output.key.fields=3",
                    "-partitioner",
                    "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
                    "-input",
                    "s3://rail-experiments/geuvadis_again_intermediate/align_reads/end_to_end_sam,s3://rail-experiments/geuvadis_again_intermediate/realign_reads/splice_sam",
                    "-output",
                    "s3://rail-experiments/geuvadis_again_intermediate/bam",
                    "-mapper",
                    "cat",
                    "-reducer",
                    "pypy /mnt/src/rail-rna/bam.py --out=S3://rail-experiments/geuvadis_again/bam --bowtie-idx=/mnt/index/genome --samtools-exe=samtools --bam-basename=alignments --manifest=/mnt/MANIFEST --keep-alive"
                ],
                "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar"
            },
            "Name": "Bam"
        }
    ]"""

"""[
{{
  "Name" : "Align reads to genome",
  "ActionOnFailure" : "CANCEL_AND_WAIT",
  "HadoopJarStep": {{
    "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar",
    "Args": [
      "-libjars", "/mnt/lib/multiplefiles.jar",
      "-D", "mapred.reduce.tasks=0",
      "-input", "/Users/anellore/rail/example/dmel_flux/preprocess",
      "-output", "/Users/anellore/rail/example/dmel_flux/intermediate/align_reads",
      "-mapper", "pypy /Users/anellore/rail/src/rail-rna/align_reads.py --bowtie-idx=/Users/anellore/Downloads/Drosophila_melanogaster_UCSC_dm3/Drosophila_melanogaster/UCSC/dm3/Sequence/BowtieIndex/genome --bowtie2-idx=/Users/anellore/Downloads/Drosophila_melanogaster_UCSC_dm3/Drosophila_melanogaster/UCSC/dm3/Sequence/Bowtie2Index/genome --bowtie2-exe=bowtie2 --exon-differentials --partition-length 5000 --manifest=/Users/anellore/rail/example/dmel_flux/dmel_flux.abs.manifest --verbose -- --local",
      "-reducer", "cat",
      "-outputformat", "edu.jhu.cs.MultipleOutputFormat"
    ]
  }}
}},
{{
  "Name" : "Aggregate duplicate read sequences to reduce realignment burden",
  "ActionOnFailure" : "CANCEL_AND_WAIT",
  "HadoopJarStep": {{
    "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar",
    "Args": [
      "-D", "mapred.reduce.tasks=16",
      "-D", "mapred.text.key.partitioner.options=-k1,1",
      "-D", "stream.num.map.output.key.fields=1",
      "-partitioner", "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-input", "/Users/anellore/rail/example/dmel_flux/intermediate/align_reads/readletize",
      "-output", "/Users/anellore/rail/example/dmel_flux/intermediate/combine_sequences",
      "-mapper", "cat",
      "-reducer", "pypy /Users/anellore/rail/src/rail-rna/sum.py --type 3 --value-count 2"
    ]
  }}
}},
{{
  "Name" : "Segment reads into readlets",
  "ActionOnFailure" : "CANCEL_AND_WAIT",
  "HadoopJarStep": {{
    "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar",
    "Args": [
      "-D", "mapred.reduce.tasks=16",
      "-D", "mapred.text.key.partitioner.options=-k1,1",
      "-D", "stream.num.map.output.key.fields=1",
      "-partitioner", "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-input", "/Users/anellore/rail/example/dmel_flux/intermediate/combine_sequences/",
      "-output", "/Users/anellore/rail/example/dmel_flux/intermediate/readletize",
      "-mapper", "cat",
      "-reducer", "pypy /Users/anellore/rail/src/rail-rna/readletize.py --max-readlet-size 23 --readlet-interval 4 --capping-multiplier 1.200000"
    ]
  }}
}},
{{
  "Name" : "Aggregate readlet sequences to reduce realignment burden",
  "ActionOnFailure" : "CANCEL_AND_WAIT",
  "HadoopJarStep": {{
    "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar",
    "Args": [
      "-D", "mapred.reduce.tasks=16",
      "-D", "mapred.text.key.partitioner.options=-k1,1",
      "-D", "stream.num.map.output.key.fields=1",
      "-partitioner", "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-input", "/Users/anellore/rail/example/dmel_flux/intermediate/readletize/",
      "-output", "/Users/anellore/rail/example/dmel_flux/intermediate/combine_subsequences",
      "-mapper", "cat",
      "-reducer", "pypy /Users/anellore/rail/src/rail-rna/sum.py --type 3"
    ]
  }}
}},
{{
  "Name" : "Align readlets to genome",
  "ActionOnFailure" : "CANCEL_AND_WAIT",
  "HadoopJarStep": {{
    "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar",
    "Args": [
      "-D", "mapred.reduce.tasks=16",
      "-D", "mapred.text.key.partitioner.options=-k1,1",
      "-D", "stream.num.map.output.key.fields=1",
      "-partitioner", "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-input", "/Users/anellore/rail/example/dmel_flux/intermediate/combine_subsequences",
      "-output", "/Users/anellore/rail/example/dmel_flux/intermediate/align_readlets",
      "-mapper", "cat",
      "-reducer", "pypy /Users/anellore/rail/src/rail-rna/align_readlets.py --bowtie-idx=/Users/anellore/Downloads/Drosophila_melanogaster_UCSC_dm3/Drosophila_melanogaster/UCSC/dm3/Sequence/BowtieIndex/genome --bowtie-exe=bowtie --verbose -- -t --sam-nohead --startverbose -v 0 -a -m 80"
    ]
  }}
}},
{{
  "Name" : "Search for introns using readlet alignments",
  "ActionOnFailure" : "CANCEL_AND_WAIT",
  "HadoopJarStep": {{
    "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar",
    "Args": [
      "-D", "mapred.reduce.tasks=16",
      "-D", "mapred.text.key.partitioner.options=-k1,1",
      "-D", "stream.num.map.output.key.fields=1",
      "-partitioner", "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-input", "/Users/anellore/rail/example/dmel_flux/intermediate/align_readlets",
      "-output", "/Users/anellore/rail/example/dmel_flux/intermediate/intron_search",
      "-mapper", "cat",
      "-reducer", "pypy /Users/anellore/rail/src/rail-rna/intron_search.py --bowtie-idx=/Users/anellore/Downloads/Drosophila_melanogaster_UCSC_dm3/Drosophila_melanogaster/UCSC/dm3/Sequence/BowtieIndex/genome --partition-length 5000 --max-intron-size 500000 --min-intron-size 10 --min-exon-size 9 --search-window-size 1000 --motif-radius 5 --verbose"
    ]
  }}
}},
{{
  "Name" : "Enumerate possible intron co-occurrences on readlets",
  "ActionOnFailure" : "CANCEL_AND_WAIT",
  "HadoopJarStep": {{
    "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar",
    "Args": [
      "-D", "mapred.reduce.tasks=16",
      "-D", "mapred.text.key.partitioner.options=-k1,2",
      "-D", "stream.num.map.output.key.fields=4",
      "-partitioner", "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-input", "/Users/anellore/rail/example/dmel_flux/intermediate/intron_search",
      "-output", "/Users/anellore/rail/example/dmel_flux/intermediate/intron_config",
      "-mapper", "cat",
      "-reducer", "pypy /Users/anellore/rail/src/rail-rna/intron_config.py --readlet-size 23 --verbose"
    ]
  }}
}},
{{
  "Name" : "Obtain exonic reference sequences appropriate for readlet alignment from intron co-occurrences",
  "ActionOnFailure" : "CANCEL_AND_WAIT",
  "HadoopJarStep": {{
    "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar",
    "Args": [
      "-D", "mapred.reduce.tasks=4",
      "-D", "mapred.text.key.partitioner.options=-k1,4",
      "-D", "stream.num.map.output.key.fields=4",
      "-partitioner", "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-input", "/Users/anellore/rail/example/dmel_flux/intermediate/intron_config",
      "-output", "/Users/anellore/rail/example/dmel_flux/intermediate/intron_fasta",
      "-mapper", "cat",
      "-reducer", "pypy /Users/anellore/rail/src/rail-rna/intron_fasta.py --verbose --bowtie-idx=/Users/anellore/Downloads/Drosophila_melanogaster_UCSC_dm3/Drosophila_melanogaster/UCSC/dm3/Sequence/BowtieIndex/genome"
    ]
  }}
}},
{{
  "Name" : "Create Bowtie 1 index from exonic reference sequences",
  "ActionOnFailure" : "CANCEL_AND_WAIT",
  "HadoopJarStep": {{
    "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar",
    "Args": [
      "-D", "mapred.reduce.tasks=1",
      "-D", "mapred.text.key.partitioner.options=-k1,1",
      "-D", "stream.num.map.output.key.fields=1",
      "-partitioner", "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-input", "/Users/anellore/rail/example/dmel_flux/intermediate/intron_fasta",
      "-output", "/Users/anellore/rail/example/dmel_flux/intermediate/intron_index",
      "-mapper", "cat",
      "-reducer", "pypy /Users/anellore/rail/src/rail-rna/intron_index.py --bowtie-build-exe=bowtie-build --bowtie-idx=/Users/anellore/Downloads/Drosophila_melanogaster_UCSC_dm3/Drosophila_melanogaster/UCSC/dm3/Sequence/BowtieIndex/genome --out=local_out/index"
    ]
  }}
}},
{{
  "Name" : "Align readlets to Bowtie 1 index of exonic sequences",
  "ActionOnFailure" : "CANCEL_AND_WAIT",
  "HadoopJarStep": {{
    "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar",
    "Args": [
      "-D", "mapred.reduce.tasks=16",
      "-D", "mapred.text.key.partitioner.options=-k1,1",
      "-D", "stream.num.map.output.key.fields=1",
      "-archives", "s3n://rail-experiments/geu4/index/intron.tar.gz#intron",
      "-partitioner", "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-input", "/Users/anellore/rail/example/dmel_flux/intermediate/combine_subsequences",
      "-output", "/Users/anellore/rail/example/dmel_flux/intermediate/realign_readlets",
      "-mapper", "cat",
      "-reducer", "pypy /Users/anellore/rail/src/rail-rna/align_readlets.py --bowtie-idx=local_out/index/intron --bowtie-exe=bowtie --verbose -- -t --sam-nohead --startverbose -v 0 -a -m 80"
    ]
  }}
}},
{{
  "Name" : "Search for intron co-occurrences on reads from readlet alignments",
  "ActionOnFailure" : "CANCEL_AND_WAIT",
  "HadoopJarStep": {{
    "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar",
    "Args": [
      "-D", "mapred.reduce.tasks=16",
      "-D", "mapred.text.key.partitioner.options=-k1,1",
      "-D", "stream.num.map.output.key.fields=1",
      "-partitioner", "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-input", "/Users/anellore/rail/example/dmel_flux/intermediate/realign_readlets,/Users/anellore/rail/example/dmel_flux/intermediate/align_readlets",
      "-output", "/Users/anellore/rail/example/dmel_flux/intermediate/cointron_search",
      "-mapper", "cat",
      "-reducer", "pypy /Users/anellore/rail/src/rail-rna/cointron_search.py --verbose"
    ]
  }}
}},
{{
  "Name" : "Obtain exonic reference sequences appropriate for read alignment from intron co-occurrences",
  "ActionOnFailure" : "CANCEL_AND_WAIT",
  "HadoopJarStep": {{
    "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar",
    "Args": [
      "-D", "mapred.reduce.tasks=16",
      "-D", "mapred.text.key.partitioner.options=-k1,2",
      "-D", "stream.num.map.output.key.fields=2",
      "-partitioner", "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-input", "/Users/anellore/rail/example/dmel_flux/intermediate/cointron_search/",
      "-output", "/Users/anellore/rail/example/dmel_flux/intermediate/cointron_fasta",
      "-mapper", "cat",
      "-reducer", "pypy /Users/anellore/rail/src/rail-rna/cointron_fasta.py --verbose --bowtie-idx=/Users/anellore/Downloads/Drosophila_melanogaster_UCSC_dm3/Drosophila_melanogaster/UCSC/dm3/Sequence/BowtieIndex/genome"
    ]
  }}
}},
{{
  "Name" : "Align reads to Bowtie 2 indexes of exome",
  "ActionOnFailure" : "CANCEL_AND_WAIT",
  "HadoopJarStep": {{
    "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar",
    "Args": [
      "-D", "mapred.reduce.tasks=16",
      "-D", "mapred.text.key.partitioner.options=-k1,1",
      "-D", "stream.num.map.output.key.fields=1",
      "-libjars", "/mnt/lib/multiplefiles.jar",
      "-partitioner", "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-input", "/Users/anellore/rail/example/dmel_flux/intermediate/align_reads/unmapped,/Users/anellore/rail/example/dmel_flux/intermediate/cointron_fasta/",
      "-output", "/Users/anellore/rail/example/dmel_flux/intermediate/realign_reads",
      "-mapper", "cat",
      "-reducer", "pypy /Users/anellore/rail/src/rail-rna/realign_reads.py --original-idx=/Users/anellore/Downloads/Drosophila_melanogaster_UCSC_dm3/Drosophila_melanogaster/UCSC/dm3/Sequence/BowtieIndex/genome --bowtie2-exe=bowtie2 --partition-length 5000 --exon-differentials --manifest=/Users/anellore/rail/example/dmel_flux/dmel_flux.abs.manifest --verbose -- --end-to-end",
      "-outputformat", "edu.jhu.cs.MultipleOutputFormat"
    ]
  }}
}},
{{
  "Name" : "Merge exon differentials at same genomic location",
  "ActionOnFailure" : "CANCEL_AND_WAIT",
  "HadoopJarStep": {{
    "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar",
    "Args": [
      "-D", "mapred.reduce.tasks=16",
      "-D", "mapred.text.key.partitioner.options=-k1,3",
      "-D", "stream.num.map.output.key.fields=3",
      "-partitioner", "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-input", "/Users/anellore/rail/example/dmel_flux/intermediate/align_reads/exon_diff,/Users/anellore/rail/example/dmel_flux/intermediate/realign_reads/exon_diff",
      "-output", "/Users/anellore/rail/example/dmel_flux/intermediate/collapse",
      "-mapper", "cat",
      "-reducer", "pypy /Users/anellore/rail/src/rail-rna/sum.py"
    ]
  }}
}},
{{
  "Name" : "Compile coverage vectors from exon differentials",
  "ActionOnFailure" : "CANCEL_AND_WAIT",
  "HadoopJarStep": {{
    "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar",
    "Args": [
      "-D", "mapred.reduce.tasks=16",
      "-D", "mapred.text.key.partitioner.options=-k1,2",
      "-D", "stream.num.map.output.key.fields=3",
      "-libjars", "/mnt/lib/multiplefiles.jar",
      "-partitioner", "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-input", "/Users/anellore/rail/example/dmel_flux/intermediate/collapse/",
      "-output", "/Users/anellore/rail/example/dmel_flux/intermediate/coverage_pre",
      "-mapper", "cat",
      "-reducer", "pypy /Users/anellore/rail/src/rail-rna/coverage_pre.py --bowtie-idx=/Users/anellore/Downloads/Drosophila_melanogaster_UCSC_dm3/Drosophila_melanogaster/UCSC/dm3/Sequence/BowtieIndex/genome --partition-stats",
      "-outputformat", "edu.jhu.cs.MultipleOutputFormat"
    ]
  }}
}},
{{
  "Name" : "Write bigbed files summarizing exome coverage by sample",
  "ActionOnFailure" : "CANCEL_AND_WAIT",
  "HadoopJarStep": {{
    "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar",
    "Args": [
      "-D", "mapred.reduce.tasks=16",
      "-D", "mapred.text.key.partitioner.options=-k1,1",
      "-D", "stream.num.map.output.key.fields=3",
      "-partitioner", "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-input", "/Users/anellore/rail/example/dmel_flux/intermediate/coverage_pre/coverage",
      "-output", "/Users/anellore/rail/example/dmel_flux/intermediate/coverage",
      "-mapper", "cat",
      "-reducer", "pypy /Users/anellore/rail/src/rail-rna/coverage.py --bowtie-idx=/Users/anellore/Downloads/Drosophila_melanogaster_UCSC_dm3/Drosophila_melanogaster/UCSC/dm3/Sequence/BowtieIndex/genome --percentile 0.750000 --out=local_out/coverage --bigbed-exe=/Users/anellore/Downloads/bedToBigBed --manifest=/Users/anellore/rail/example/dmel_flux/dmel_flux.abs.manifest --verbose"
    ]
  }}
}},
{{
  "Name" : "Write normalization factors corresponding to samples",
  "ActionOnFailure" : "CANCEL_AND_WAIT",
  "HadoopJarStep": {{
    "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar",
    "Args": [
      "-D", "mapred.reduce.tasks=1",
      "-D", "mapred.text.key.partitioner.options=-k1,1",
      "-D", "stream.num.map.output.key.fields=2",
      "-partitioner", "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-input", "/Users/anellore/rail/example/dmel_flux/intermediate/coverage",
      "-output", "/Users/anellore/rail/example/dmel_flux/intermediate/coverage_post",
      "-mapper", "cat",
      "-reducer", "pypy /Users/anellore/rail/src/rail-rna/coverage_post.py --out=local_out/normalize --manifest=/Users/anellore/rail/example/dmel_flux/dmel_flux.abs.manifest"
    ]
  }}
}},
{{
  "Name" : "Aggregate intron and indel results for each sample",
  "ActionOnFailure" : "CANCEL_AND_WAIT",
  "HadoopJarStep": {{
    "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar",
    "Args": [
      "-D", "mapred.reduce.tasks=16",
      "-D", "mapred.text.key.partitioner.options=-k1,6",
      "-D", "stream.num.map.output.key.fields=6",
      "-partitioner", "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-input", "/Users/anellore/rail/example/dmel_flux/intermediate/realign_reads/bed,/Users/anellore/rail/example/dmel_flux/intermediate/align_reads/bed",
      "-output", "/Users/anellore/rail/example/dmel_flux/intermediate/bed_pre",
      "-mapper", "cat",
      "-reducer", "pypy /Users/anellore/rail/src/rail-rna/bed_pre.py"
    ]
  }}
}},
{{
  "Name" : "Write bed files with intron and indel results by sample",
  "ActionOnFailure" : "CANCEL_AND_WAIT",
  "HadoopJarStep": {{
    "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar",
    "Args": [
      "-D", "mapred.reduce.tasks=1",
      "-D", "mapred.text.key.partitioner.options=-k1,2",
      "-D", "stream.num.map.output.key.fields=4",
      "-partitioner", "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-input", "/Users/anellore/rail/example/dmel_flux/intermediate/bed_pre/",
      "-output", "/Users/anellore/rail/example/dmel_flux/intermediate/bed",
      "-mapper", "cat",
      "-reducer", "pypy /Users/anellore/rail/src/rail-rna/bed.py --bowtie-idx=/Users/anellore/Downloads/Drosophila_melanogaster_UCSC_dm3/Drosophila_melanogaster/UCSC/dm3/Sequence/BowtieIndex/genome  --out=local_out/bed --manifest=/Users/anellore/rail/example/dmel_flux/dmel_flux.abs.manifest --bed-basename="
    ]
  }}
}},
{{
  "Name" : "Write bam files containing final read alignments by sample",
  "ActionOnFailure" : "CANCEL_AND_WAIT",
  "HadoopJarStep": {{
    "Jar": "/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar",
    "Args": [
      "-D", "mapred.reduce.tasks=1",
      "-D", "mapred.text.key.partitioner.options=-k1,1",
      "-D", "stream.num.map.output.key.fields=3",
      "-partitioner", "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-input", "/Users/anellore/rail/example/dmel_flux/intermediate/end_to_end_sam,/Users/anellore/rail/example/dmel_flux/intermediate/realign_reads/splice_sam",
      "-output", "/Users/anellore/rail/example/dmel_flux/intermediate/bam",
      "-mapper", "cat",
      "-reducer", "pypy /Users/anellore/rail/src/rail-rna/bam.py --out=local_out/bam --bowtie-idx=/Users/anellore/Downloads/Drosophila_melanogaster_UCSC_dm3/Drosophila_melanogaster/UCSC/dm3/Sequence/BowtieIndex/genome  --samtools-exe=samtools --bam-basename=alignments --manifest=/Users/anellore/rail/example/dmel_flux/dmel_flux.abs.manifest"
    ]
  }}
}}
]""".format(

	)

    def preprocess_emr_json(self, json_file):

    def process_emr_json(self, json_file):

    def preprocess_step_json(self, json_file):

    def process_step_json(self, json_file):


