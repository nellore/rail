"""
rail-rna_config.py
Part of Rail-RNA

Contains classes that perform error checking and generate JSON configuration
output for Rail-RNA. These configurations are parsable by Dooplicity's
hadoop_simulator.py and emr_runner.py.

Class structure is designed so only those arguments relevant to modes/job flows
are included.

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
utils_path = os.path.join(base_path, 'rna', 'utils')
site.addsitedir(utils_path)
site.addsitedir(base_path)
import dooplicity.ansibles as ab
import tempfile
import shutil
from tools import which, is_exe, path_join
from argparse import SUPPRESS

def step(name, inputs, output, mapper='cat', reducer='cat', 
    action_on_failure='TERMINATE_JOB_FLOW',
    jar='/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar',
    tasks=0, partitioner_options=None, key_fields=None, archives=None,
    multiple_outputs=False, inputformat=None):
    """ Outputs JSON for a given step.

        name: name of step
        inputs: list of input directories/files
        output: output directory
        mapper: mapper command
        reducer: reducer command
        jar: path to Hadoop Streaming jar; ignored in local mode
        tasks: reduce task count
        partitioner options: UNIX sort-like partitioner options
        key fields: number of key fields,
        archives: -archives option
        multiple_outputs: True iff there are multiple outputs; else False
        inputformat: -inputformat option

        Return value: step dictionary
    """
    to_return = {
        'Name' : name
        'ActionOnFailure' : action_on_failure,
        'HadoopJarStep' : {
            'Jar' : jar,
            'Args' : []
        }

    }
    if multiple_outputs:
        # This only matters on EMR
        to_return['HadoopJarStep']['Args'].extend([
            '-libjars', '/mnt/lib/multiplefiles.jar'
        ])
    to_return.extend(['-D', 'mapred.reduce.tasks=%d' % tasks])
    if partioner_options is not None and key_fields is not None:
        to_return['HadoopJarStep']['Args'].extend([
                '-D', 'mapred.text.key.partitioner.options=-%s'
                            % partitioner_options,
                '-D', 'stream.num.map.output.key.fields=%d' % key_fields
            ])
    if archives is not None:
        to_return['HadoopJarStep']['Args'].extend([
                '-archives', archives
            ])
    to_return['HadoopJarStep']['Args'].extend([
            '-partitioner',
            'org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner',
        ])
    for an_input in inputs.split(',')
        to_return['HadoopJarStep']['Args'].extend([
                '-input', an_input.strip()
            ])
    to_return['HadoopJarStep']['Args'].extend([
            '-output', output,
            '-mapper', mapper,
            '-reducer', reducer
        ])
    if multiple_outputs:
        to_return.extend([
                '-outputformat', 'edu.jhu.cs.MultipleOutputFormat'
            ])
    if input_format is not None:
        to_return.extend([
                '-intputformat', inputformat
            ])
    return to_return

def steps(pseudosteps, action_on_failure, jar, step_dir, 
            reducer_count, intermediate_dir, unix=False):
    """ Turns list with "pseudosteps" into well-formed StepConfig list.

        A pseudostep looks like this:

            {
                'name' : [name of step]
                'run' : Python script name; like 'preprocess.py' + args
                'inputs' : list of input directories
                'no_input_prefix' : key that's present iff intermediate dir
                    should not be prepended to inputs
                'output' : output directory
                'no_output_prefix' : key that's present iff intermediate dir
                    should not be prepended to output dir
                'keys'  : Number of key fields; present only if reducer
                'part'  : KeyFieldBasedPartitioner options; present only if
                            reducer
                'taskx' : number of tasks per reducer
                'inputformat' : input format; present only if necessary
                'archives' : archives parameter; present only if necessary
                'multiple_outputs' : key that's present iff there are multiple
                    outputs
            }

        pseudosteps: array of pseudosteps
        action_on_failure: action on failure to take
        jar: path to Hadoop Streaming jar
        step_dir: where to find Python scripts for steps
        reducer_count: number of reducers; determines number of tasks
        unix: performs UNIX-like path joins

        Return value: list of StepConfigs (see Elastic MapReduce API docs)
    """
    true_steps = []
    for pseudostep in pseudosteps:
        assert ('keys' in pseudostep and 'part' in pseudostep) or
                ('keys' not in pseudostep and 'part' not in pseudostep)
        true_steps.append(step(
                            name=pseudostep['name'],
                            inputs=([path_join(unix, intermediate_dir,
                                        an_input) for an_input in
                                        pseudostep['inputs']]
                                    if 'no_input_prefix' not in
                                    pseudostep else pseudostep['inputs']),
                            output=(path_join(unix, intermediate_dir,
                                                    pseudostep['output'])
                                    if 'no_output_prefix' not in
                                    pseudostep else pseudostep['output']),
                            mapper=(path_join(unix, step_dir,
                                                    pseudostep['run'])
                                    if 'keys' not in pseudostep
                                    else 'cat'),
                            reducer=(path_join(unix, step_dir,
                                                    pseudostep['run'])
                                    if 'keys' in pseudostep
                                    else 'cat')
                            action_on_failure=action_on_failure,
                            jar=jar,
                            tasks=reducer_count * pseudostep['taskx'],
                            partitioner_options=(pseudostep['part']
                                if 'part' in pseudostep else None),
                            key_fields=(pseudostep['keys']
                                if 'keys' in pseudostep else None),
                            archives=(pseudostep['archives']
                                if 'archives' in pseudostep else None)
                            multiple_outputs=(True if 'multiple_outputs'
                                    in pseudostep else False
                                )
                            inputformat=(pseudostep['inputformat']
                                if 'inputformat' in pseudostep else None)
                        )
                    )
    return true_steps

class RailRnaErrors:
    """ Holds accumulated errors in Rail-RNA's input parameters.

        Checks only those parameters common to all modes/job flows.
    """
    def __init__(self, manifest, output_dir,
            intermediate_dir='./intermediate', force=False, aws_exe=None,
            profile='default', region='us-east-1', verbose=False
        ):
        '''Store all errors uncovered in a list, then output. This prevents the
        user from having to rerun Rail-RNA to find what else is wrong with
        the command-line parameters.'''
        self.errors = []
        self.manifest_dir = None
        self.manifest = manifest
        self.output_dir = output_dir
        self.intermediate_dir = intermediate_dir
        self.aws_exe = aws_exe
        self.region = region
        self.force = force
        self.checked_programs = set()
        self.curl_exe = curl_exe
        self.verbose = verbose

    def check_s3(self, reason=None):
        """ Checks for AWS CLI and configuration file.

            In this script, S3 checking is performed as soon as it is found
            that S3 is needed. If anything is awry, a RuntimeError is raised
            _immediately_ (the standard behavior is to raise a RuntimeError
            only after errors are accumulated). A reason specifying where
            S3 credentials were first needed can also be provided.

            reason: string specifying where S3 credentials were first
                needed.

            No return value.
        """
        original_errors_size = len(self.errors)
        if aws_exe is None:
            self.aws_exe = 'aws'
            if not which(self.aws_exe):
                self.errors.append(('The AWS CLI executable '
                                    'was not found. Make sure that the '
                                    'executable is in PATH, or specify the '
                                    'location of the executable with '
                                    '"--aws-exe".'))
            else:
                self.errors.append(('The AWS CLI executable ("--aws-exe") '
                                    '"{0}" was not found. Make sure that '
                                    'the file is present and is '
                                    'executable.').format(aws_exe))
        elif not is_exe(self.aws_exe):
            self.errors.append(('The AWS CLI executable ("--aws-exe") '
                                '"{0}" is either not present or not '
                                'executable.').format(aws_exe))
        self._aws_access_key_id = None
        self._aws_secret_access_key = None
        if profile == 'default':
            # Search environment variables for keys first if profile is default
            try:
                self._aws_access_key_id = os.environ['AWS_ACCESS_KEY_ID']
                self._aws_secret_access_key \
                    = os.environ['AWS_SECRET_ACCESS_KEY']
                to_search = None
            except KeyError:
                to_search = '[default]'
            try:
                # Also grab region
                self.region = os.environ['AWS_DEFAULT_REGION']
            except KeyError:
                pass
        else:
            to_search = '[profile ' + profile + ']'
        # Now search AWS CLI config file for the right profile
        if to_search is not None:
            config_file = os.path.join(os.environ['HOME'], '.aws', 'config')
            try:
                with open(config_file) as config_stream:
                    for line in config_stream:
                        if line.strip() == to_search:
                            break
                    for line in config_stream:
                        tokens = [token.strip() for token in line.split('=')]
                        if tokens[0] == 'region' \
                            and self.region == 'us-east-1':
                            self.region = tokens[1]
                        elif tokens[0] == 'aws_access_key_id':
                            self._aws_access_key_id = tokens[1]
                        elif tokens[0] == 'aws_secret_access_key':
                            self._aws_secret_access_key = tokens[1]
                        else:
                            line = line.strip()
                            if line[0] == '[' and line[-1] == ']':
                                # Break on start of new profile
                                break
            except IOError:
                self.errors.append(
                                   ('No valid AWS CLI configuration found. '
                                    'Make sure the AWS CLI is installed '
                                    'properly and that one of the following '
                                    'is true:\n\na) The environment variables '
                                    '"AWS_ACCESS_KEY_ID" and '
                                    '"AWS_SECRET_ACCESS_KEY" are set to '
                                    'the desired AWS access key ID and '
                                    'secret access key, respectively, and '
                                    'the profile ("--profile") is set to '
                                    '"default" (its default value).\n\n'
                                    'b) The file ".aws/config" exists in your '
                                    'home directory with a valid profile. '
                                    'To set this file up, run "aws --config" '
                                    'after installing the AWS CLI.')
                                )
        if len(self.errors) != original_errors_size:
            if reason:
                raise RuntimeError(('\n'.join(['%d) %s' % (i, error)
                                    for i, error
                                    in enumerate(self.errors)]) + 
                                    '\n\nNote that the AWS CLI is needed '
                                    'because {0}. If all dependence on S3 in '
                                    'the pipeline is removed, the AWS CLI '
                                    'need not be installed.').format(reason))
            else:
                raise RuntimeError(('\n'.join(['%d) %s' % (i, error)
                                    for i, error
                                    in enumerate(self.errors)]) + 
                                    '\n\nIf all dependence on S3 in the '
                                    'pipeline is removed, the AWS CLI need '
                                    'not be installed.'))
        self.checked_programs.add('AWS CLI')

    def check_program(exe, program_name, parameter, var_to_set,
                        entered_exe=None, reason=None):
        """ Checks if program in PATH or if user specified it properly.

            Errors are added to self.errors.

            exe: executable to search for
            program name: name of program
            parameter: corresponding command line parameter
                (e.g., --bowtie-exe)
            var_to_set: variable to set (e.g., "object.bowtie_exe")
            entered_exe: None if the user didn't enter an executable; otherwise
                whatever the user entered
            reason: FOR CURL ONLY: raise RuntimeError _immediately_ if Curl
                not found but needed

            No return value.
        """
        original_errors_size = len(self.errors)
        if entered_exe is None:
            if not which(exe):
                self.errors.append(
                        ('The executable "{0}" for {1} was either not found '
                         'in PATH or is not executable. Check that the '
                         'program is installed properly and executable; then '
                         'either add the executable to PATH or specify it '
                         'directly with "{3}".').format(exe, program_name,
                                                            parameter)
                    )
            else:
                var_to_set = exe
        elif not is_exe(entered_exe):
            self.errors.append(
                    ('The executable "{0}" entered for {1} via "{2}" was '
                     'either not found or is not executable.').format(exe,
                                                                program_name,
                                                                parameter)
                )
        else:
            var_to_set = entered_exe
        if original_errors_size != len(self.errors) and reason:
            raise RuntimeError(('\n'.join(['%d) %s' % (i, error)
                                for i, error
                                in enumerate(self.errors)]) + 
                                '\n\nNote that Curl is needed because {0}.'
                                ' If all dependence on web resources is '
                                'removed from the pipeline, Curl need '
                                'not be installed.').format(reason))
        self.checked_programs.add(program_name)

    @staticmethod
    def add_args(parser):
        parser.add_argument(
            '--aws-exe', type=str, required=False,
            default=None,
            help=('AWS CLI executable. If "aws" is in PATH, this parameter '
                 'should be left unspecified. If it\'s not, the full path to '
                 'the executable should be specified. If S3 is never used in '
                 'a given job flow, this parameter is inconsequential.')
        )
        parser.add_argument(
            '--curl-exe', type=str, required=False,
            default=None,
            help=('Curl executable. If "curl" is in PATH, this parameter '
                  'should be left unspecified. If it\'s not, the full path to '
                  'the executable should be specified.')
        )
        parser.add_argument(
            '--profile', type=str, required=False,
            default=None,
            help=('AWS CLI profile to use. Defaults to [default]; if, '
                  'however, the environment variables "AWS_ACCESS_KEY_ID", '
                  'and "AWS_SECRET_ACCESS_KEY" are set, these are used '
                  'instead of [default].')
        )
        parser.add_argument(
            '-f', '--force', action='store_const', const=True,
            default=False,
            help='Overwrites output directory if it exists.'
        )
        '''--region's help looks different from mode to mode; don't include it
        here.'''

class RailRnaLocal:
    """ Checks local-mode JSON from input parameters and relevant programs.

        Subsumes only those parameters relevant to local mode. Adds errors
        to base instance of RailRnaErrors.
    """
    def __init__(self, base, check_manifest=False,
                    num_processes=1, keep_intermediates=False):
        """ base: instance of RailRnaErrors """
        # Initialize ansible for easy checks
        ansible = ab.Ansible()             )
        if not ab.Url(base,intermediate_dir).is_local:
            base.errors.append(('Intermediate directory must be local '
                                'when running Rail-RNA in local ("--local") '
                                'mode, but {0} was entered.').format(
                                        base.intermediate_dir
                                    ))
        output_dir_url = ab.Url(base.output_dir)
        if output_dir_url.is_curlable:
            base.errors.append(('Output directory must be local or on S3 '
                                'when running Rail-RNA in local ("--local") '
                                'mode, but {0} was entered.').format(
                                        base.output_dir
                                    ))
        elif output_dir_url.is_s3 and 'AWS CLI' not in base.checked_programs:
            base.check_s3(reason='the output directory is on S3')
            # Change ansible params
            ansible.aws_exe = base.aws_exe
            ansible.profile = base.profile
        if not base.force:
            if output_dir_url.is_local \
                and os.path.exists(output_dir_url.to_url()):
                base.errors.append(('Output directory {0} exists, '
                                    'and "--force" was not invoked to permit '
                                    'overwriting it.').format(base.output_dir))
            elif output_dir_url.is_s3 \
                and ansible.s3_ansible.is_dir(base.output_dir):
                base.errors.append(('Output directory {0} exists on S3, and '
                                    '"--force" was not invoked to permit '
                                    'overwriting it.').format(base_output_dir))
        # Check manifest; download it if necessary
        manifest_url = ab.Url(base.manifest)
        if manifest_url.is_s3 and 'AWS CLI' not in base.checked_programs:
            base.check_s3(reason='the manifest file is on S3')
            # Change ansible params
            ansible.aws_exe = base.aws_exe
            ansible.profile = base.profile
        elif manifest_url.is_curlable \
             and 'Curl' not in base.checked_programs:
             base.check_program('curl', 'Curl', '--curl_exe', base.curl_exe,
                                    entered_exe=base.curl_exe,
                                    reason='the manifest file is on the web')
             ansible.curl_exe = base.curl_exe
        if not ansible.exists(manifest_url.to_url()):
            base.errors.append(('Manifest file ("--manifest") {0} '
                                'does not exist. Check the URL and '
                                'try again.').format(base.manifest))
        else:
            if not manifest_url.is_local:
                base.manifest_dir = tempfile.mkdtemp()
                base.manifest = os.path.join(base.manifest_dir, 'MANIFEST')
                ansible.get(manifest_url, destination=base.manifest)
            else:
                base.manifest = manifest
            files_to_check = []
            with open(base.manifest) as manifest_stream:
                for line in manifest_stream:
                    tokens = line.strip().split('\t')
                    if len(tokens) == 5:
                        files_to_check.extend([tokens[0], tokens[2]])
                    elif len(tokens) == 3:
                        files_to_check.append(tokens[0])
                    else:
                        base.errors.append(('The following line from the '
                                            'manifest file {0} '
                                            'has an invalid number of '
                                            'tokens:\n{1}'
                                            ).format(
                                                    manifest,
                                                    line
                                                ))
            if files_to_check and check_manifest:
                # Check files in manifest only if in preprocess job flow
                for filename in files_to_check:
                    filename_url = ab.Url(filename)
                    if filename_url.is_s3 \
                        and 'AWS CLI' not in base.checked_programs:
                            base.check_s3(reason=('at least one sample '
                                                  'FASTA/FASTQ from the '
                                                  'manifest file is on S3'))
                            # Change ansible params
                            ansible.aws_exe = base.aws_exe
                            ansible.profile = base.profile
                    elif filename_url.is_curlable \
                        and 'Curl' not in base.checked_programs:
                            base.check_program('curl', 'Curl', '--curl_exe', 
                                                base.curl_exe,
                                                entered_exe=base.curl_exe,
                                                reason=('at least one sample '
                                                  'FASTA/FASTQ from the '
                                                  'manifest file is on '
                                                  'the web'))
                            ansible.curl_exe = base.curl_exe
                    if not ansible.exists(filename_url):
                        base.errors.append(('The file {0} from the manifest '
                                            'file {1} does not exist. Check '
                                            'the URL and try again.').format(
                                                                filename,
                                                                manifest
                                                            ))
            else:
                base.errors.append(('Manifest file ("--manifest") {0} '
                                    'has no valid lines.').format(
                                                                manifest
                                                            ))
        from multiprocessing import cpu_count
        if num_processes:
            if not (isinstance(num_processes, int)
                                    and num_processes >= 1):
                base.errors.append('Number of processes ("--num-processes") '
                                   'must be an integer >= 1, '
                                   'but {0} was entered.'.format(
                                                    num_processes
                                                ))
            else:
                base.num_processes = num_processes
        else:
            try:
                base.num_processes = cpu_count()
            except NotImplementedError:
                base.num_processes = 1
            if base.num_processes != 1:
                '''Make default number of processes cpu count less 1
                so Facebook tab in user's browser won't go all unresponsive.'''
                base.num_processes -= 1
        self.keep_intermediates = keep_intermediates

    @staticmethod
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
            help=('Output directory, which in local ("--local") mode must be '
                  'on the local filesystem or on S3. This directory is not '
                  'overwritten unless the force overwrite ("--force") '
                  'parameter is also invoked.')
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
            '--num-processes', type=int, required=False,
            default=None,
            help=('Number of processes to run simultaneously. This defaults '
                  'to the number of cores on the machine less 1 if more than '
                  'one core is available, or simply 1 if the program could '
                  'not determine the number of available cores.')
        )
        parser.add_argument(
            '--keep-intermediates', action='store_const', const=True,
            default=False,
            help='Keeps intermediate files after a job flow is completed.'
        )
        parser.add_argument(
            '--verbose', action='store_const', const=True,
            default=False,
            help='Outputs extra debugging statements to stderr.'
        )

class RailRnaCloud:
    """ Checks cloud-mode input parameters and relevant programs.

        Subsumes only those parameters relevant to cloud mode. Adds errors
        to base instance of RailRnaErrors.
    """
    def __init__(self, base, log_uri=None, ami_version='2.4.2',
        visible_to_all_users=False, tags='',
        name='Rail-RNA Job Flow',
        action_on_failure='TERMINATE_JOB_FLOW',
        hadoop_jar='/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar',
        master_instance_count=1, master_instance_type='c1.xlarge',
        master_instance_bid_price=None, core_instance_count=1,
        core_instance_type=None, core_instance_bid_price=None,
        task_instance_count=0, task_instance_type=None,
        task_instance_bid_price=None, ec2_key_name=None, keep_alive=False,
        termination_protected=False):

        # CLI is REQUIRED in cloud mode
        base.check_s3(reason='Rail-RNA is running in cloud ("--cloud") mode')

        # Initialize possible options
        instance_core_counts = {
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

        instance_swap_allocations = {
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
        32- vs. 64-bit issues: instance_bits = {
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

        if log_uri is not None and not Url(log_uri).is_s3:
            base.errors.append('Log URI ("--log-uri") must be on S3, but '
                               '"{0}" was entered.'.format(log_uri))
        base.log_uri = log_uri
        base.visible_to_all_users = visible_to_all_users
        base.tags = str([tag.strip() for tag in tags.split(',')])
        base.name = name

        actions_on_failure \
            = set(['TERMINATE_JOB_FLOW', 'CANCEL_AND_WAIT', 'CONTINUE',
                    'TERMINATE_CLUSTER'])

        if action_on_failure not in actions_on_failure:
            base.errors.append('Action on failure ("--action-on-failure") '
                               'must be one of {"TERMINATE_JOB_FLOW", '
                               '"CANCEL_AND_WAIT", "CONTINUE", '
                               '"TERMINATE_CLUSTER"}, but '
                               '{0} was entered.'.format(
                                                action_on_failure
                                            ))
        base.action_on_failure = action_on_failure
        base.hadoop_jar = hadoop_jar
        if not (isinstance(num_processes, int)
                and num_processes >= 1):
            base.errors.append('Number of processes ("--num-processes") must '
                               'be an integer >= 1, '
                               'but {0} was entered.'.format(
                                                num_processes
                                            ))
        base.tasks_per_reducer = tasks_per_reducer
        base.reducer_count = reducer_count
        instance_type_message = ('Instance type ("--instance-type") must be '
                                 'in the set {"m1.small", "m1.large", '
                                 '"m1.xlarge", "c1.medium", "c1.xlarge", '
                                 '"m2.xlarge", "m2.2xlarge", "m2.4xlarge", '
                                 '"cc1.4xlarge"}, but {0} was entered.')
        if master_instance_type not in instance_core_counts:
            base.errors.append(('Master instance type '
                               '("--master-instance-type") not valid. %s')
                                % instance_type_message.format(
                                                        master_instance_type
                                                    ))
        base.master_instance_type = master_instance_type
        if core_instance_type is None:
            base.core_instance_type = base.master_instance_type
        else:
            if core_instance_type not in instance_core_counts:
                base.errors.append(('Core instance type '
                                    '("--core-instance-type") not valid. %s')
                                    % instance_type_message.format(
                                                        core_instance_type
                                                    ))
            base.core_instance_type = core_instance_type
        if task_instance_type is None:
            base.task_instance_type = base.master_instance_type
        else:
            if task_instance_type not in instance_core_counts:
                base.errors.append(('Task instance type '
                                    '("--task-instance-type") not valid. %s')
                                    % instance_type_message.format(
                                                        task_instance_type
                                                    ))
            base.task_instance_type = task_instance_type
        if master_instance_bid_price is None:
            base.spot_master = False
        else:
            if not (isinstance(master_instance bid_price, float) 
                    and master_instance_bid_price > 0):
                base.errors.append('Spot instance bid price for master nodes '
                                   '(--master-instance-bid-price) must be '
                                   '> 0, but {0} was entered.'.format(
                                                    master_instance_bid_price
                                                ))
            base.spot_master = True
            base.master_instance_bid_price = master_instance_bid_price
        if core_instance_bid_price is None:
            base.spot_core = False
        else:
            if not (isinstance(core_instance bid_price, float) 
                    and core_instance_bid_price > 0):
                base.errors.append('Spot instance bid price for core nodes '
                                   '(--core-instance-bid-price) must be '
                                   '> 0, but {0} was entered.'.format(
                                                    core_instance_bid_price
                                                ))
            base.spot_core = True
            base.core_instance_bid_price = core_instance_bid_price
        if task_instance_bid_price is None:
            base.spot_task = False
        else:
            if not (isinstance(task_instance bid_price, float) 
                    and task_instance_bid_price > 0):
                base.errors.append('Spot instance bid price for task nodes '
                                   '(--task-instance-bid-price) must be '
                                   '> 0, but {0} was entered.'.format(
                                                    task_instance_bid_price
                                                ))
            base.spot_task = True
            base.task_instance_bid_price = task_instance_bid_price
        if not (isinstance(master_instance_count, int)
                and master_instance_count >= 1):
            base.errors.append('Master instance count '
                               '("--master-instance-count") must be an '
                               'integer >= 1, but {0} was entered.'.format(
                                                    master_instance_count
                                                ))
        base.master_instance_count = master_instance_count
        if not (isinstance(core_instance_count, int)
                and core_instance_count >= 0):
            base.errors.append('Core instance count '
                               '("--core-instance-count") must be an '
                               'integer >= 1, but {0} was entered.'.format(
                                                    core_instance_count
                                                ))
        base.core_instance_count = core_instance_count
        if not (isinstance(task_instance_count, int)
                and task_instance_count >= 0):
            base.errors.append('Task instance count '
                               '("--task-instance-count") must be an '
                               'integer >= 1, but {0} was entered.'.format(
                                                    task_instance_count
                                                ))
        base.task_instance_count = task_instance_count
        if base.core_instance_count > 0:
            base.swap_allocation \
                = base.instance_swap_allocations[base.core_instance_type]
        else:
            base.swap_allocation \
                = base.instance_swap_allocations[base.master_instance_type]
        base.ec2_key_name = ec2_key_name
        base.keep_alive = keep_alive
        base.termination_protected = termination_protected

    @staticmethod
    def add_args(parser, usage=None):
        """ usage: argparse.SUPPRESS if advanced options should be suppressed;
                else None
        """
        parser.add_argument('--log-uri', type=str, required=False,
            default=None,
            usage=usage,
            help=('Directory on S3 in which to store Hadoop logs. Defaults '
                  'to "logs" subdirectory of output directory.')
        )
        parser.add_argument('--ami-version', type=str, required=False,
            default='2.4.2',
            usage=usage,
            help='Version of Amazon Linux AMI to use on EC2.'
        )
        parser.add_argument('--visible-to-all-users', action='store_const',
            const=True,
            default=False,
            usage=usage,
            help='Makes EC2 cluster visible to all IAM users within the ' \
                 'EMR CLI'
        )
        parser.add_argument('--action-on-failure', type=str, required=False,
            default='TERMINATE_JOB_FLOW',
            usage=usage,
            help='Specifies what action to take if a job flow fails on a ' \
                 'given step. Options are {"TERMINATE_JOB_FLOW", ' \
                 '"CANCEL_AND_WAIT", "CONTINUE", "TERMINATE_CLUSTER"}.'
        )
        parser.add_argument('--hadoop-jar', type=str, required=False,
            default='/home/hadoop/contrib/streaming/hadoop-streaming-1.0.3.jar'
            usage=usage,
            help='Hadoop Streaming Java ARchive to use. Controls version ' \
                 'of Hadoop Streaming.'
        )
        parser.add_argument('--master-instance-count', type=str,
            required=False,
            default=1,
            usage=usage,
            help=('Number of master instances. A master instance helps manage '
                  'the cluster, tracking the status of each task.')
        )
        parser.add_argument('-c', '--core-instance-count', type=str,
            required=False,
            default=1,
            usage=usage,
            help=('Number of core instances. A core instance runs Hadoop '
                  'maps and reduces and stores intermediate data.')
        )
        parser.add_argument('--task-instance-count', type=str,
            required=False,
            default=0,
            help=('Number of task instances. A task instance runs Hadoop '
                  'maps and reduces and but does not store any data. This is '
                  'useful if running task instances as spot instances; if '
                  'the user loses spot instances because her bid price '
                  'fell below market value, her job flow will not fail.')
        )
        parser.add_argument('--master-instance-bid-price', type=str,
            required=False,
            default=None,
            usage=usage,
            help=('Bid price for master instances (in dollars/hour). Invoke '
                  'only if running master instances as spot instances.')
        )
        parser.add_argument('--core-instance-bid-price', type=str,
            required=False,
            default=None,
            help=('Bid price for core instances (in dollars/hour). Invoke '
                  'only if running core instances as spot instances.')
        )
        parser.add_argument('--task-instance-bid-price', type=str,
            required=False,
            default=None,
            help=('Bid price for each task instance (in dollars/hour). Invoke '
                  'only if running task instances as spot instances.')
        )
        parser.add_argument('--master-instance-type', type=str,
            required=False,
            default='c1.xlarge',
            usage=usage,
            help=('Master instance type. c1.xlarge is most cost-effective '
                  'across the board for Rail-RNA.')
        )
        parser.add_argument('--core-instance-type', type=str,
            required=False,
            default=None,
            usage=usage,
            help=('Core instance type. c1.xlarge seems most cost-effective '
                  'across the board for Rail-RNA. Defaults to master '
                  'instance type if left unspecified.')
        )
        parser.add_argument('--task-instance-type', type=str,
            required=False,
            default=None,
            usage=usage,
            help=('Task instance type. c1.xlarge seems most cost-effective '
                  'across the board for Rail-RNA. Defaults to master '
                  'instance type if left unspecified.')
        )
        parser.add_argument('--ec2-key-name', type=str,
            reqired=False,
            default=None,
            usage=usage,
            help=('Name of key pair for connecting to EC2 instances via, '
                  'for example, SSH. May be useful for debugging.')
        )
        parser.add_argument('--keep-alive', type=str,
            required=False,
            default=False,
            usage=usage,
            help='Keeps EC2 cluster alive after job flow is completed.'
        )
        parser.add_argument('--termination-protected', type=str,
            required=False,
            default=False,
            usage=usage,
            help='Protects cluster from termination in case of job flow ' \
                 'failure.'
        )

class RailRnaPreprocess:
    """ Sets parameters relevant to just the preprocessing step of a job flow.
    """
    def __init__(self, base, nucleotides_per_input=8000000, gzip_input=True):
        if not (isinstance(nucleotides_per_input, int) and
                nucleotides_per_input > 0):
            base.errors.append('Nucleotides per input '
                               '(--nucleotides-per-input) must be an integer '
                               '> 0, but {0} was entered.'.format(
                                                        nucleotides_per_input
                                                       ))
        base.nucleotides_per_input = nucleotides_per_input
        base.gzip_input = gzip_input

    @staticmethod
    def add_args(parser):
        """ Adds parameter descriptions relevant to preprocess job flow to an
            object of class argparse.ArgumentParser.

            No return value.
        """
        parser.add_argument(
            '--nucleotides-per-input', type=int, required=False,
            default=8000000,
            help=('Maximum number of nucleotides from a given sample to '
                  'assign to each task. Keep this value small enough that '
                  'there are at least as many tasks as there are processor '
                  'cores available.')
        )
        parser.add_argument(
            '--do-not-gzip-output', action='store_const', const=True,
            default=False,
            help=('Leaves output of preprocess step uncompressed. This makes '
                  'preprocessing faster but takes up more hard drive space.')
        )

class RailRnaAlign:
    """ Sets parameters relevant to just the "align" job flow. """
    def __init__(self, base, local=False, bowtie1_exe=None,
        bowtie1_idx='genome', bowtie1_build_exe=None, bowtie2_exe=None,
        bowtie2_build_exe=None, bowtie2_idx='genome',
        bowtie2_args='', samtools_exe=None, bedtobigbed_exe=None,
        genome_partition_length=5000, max_readlet_size=25,
        min_readlet_size=15, readlet_interval=4, cap_size_multiplier=1.2,
        max_intron_size=500000, min_intron_size=10, min_exon_size=9,
        motif_search_window_size=1000, motif_radius=5,
        normalize_percentile=0.75, do_not_output_bam_by_chr=False,
        output_sam=False, bam_basename='alignments', bed_basename='',
        assembly='hg19', s3_ansible=None):
        if local:
            '''Programs and Bowtie indices should be checked only in local
            mode.'''
            base.check_program('bowtie', 'Bowtie 1', '--bowtie1-exe',
                                base.bowtie_exe, entered_exe=bowtie1_exe)
            base.check_program('bowtie-build', 'Bowtie 1 Build',
                                '--bowtie1-build-exe',
                                base.bowtie1_build_exe,
                                entered_exe=bowtie1_build_exe)
            for extension in ['.1.ebwt', '.2.ebwt', '.3.ebwt', '.4.ebwt', 
                                '.rev.1.ebwt', '.rev.2.ebwt']:
                index_file = bowtie1_idx + extension
                if not ab.Url(index_file).is_local:
                    base_errors.append(('Bowtie 1 index file {0} must be '
                                        'on the local filesystem.').format(
                                            index_file
                                        ))
                elif not os.path.exists(index_file):
                    base.errors.append(('Bowtie 1 index file {0} does not '
                                        'exist.').format(index_file))
            base.check_program('bowtie2', 'Bowtie 2', '--bowtie2-exe',
                                base.bowtie2_exe, entered_exe=bowtie2_exe)
            base.check_program('bowtie2-build', 'Bowtie 2 Build',
                                '--bowtie2-build-exe',
                                base.bowtie2_build_exe,
                                entered_exe=bowtie2_build_exe)
            for extension in ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', 
                                '.rev.1.bt2', '.rev.2.bt2']:
                index_file = bowtie2_idx + extension
                if not ab.Url(index_file).is_local:
                    base_errors.append(('Bowtie 2 index file {0} must be '
                                        'on the local filesystem.').format(
                                            index_file
                                        ))
                elif not os.path.exists(index_file):
                    base.errors.append(('Bowtie 2 index file {0} does not '
                                        'exist.').format(index_file))
            base.check_program('samtools', 'SAMTools', '--samtools-exe',
                                base.samtools_exe, entered_exe=samtools_exe)
            base.check_program('bedToBigBed', 'BedToBigBed',
                                '--bedtobigbed-exe', base.bedtobigbed_exe,
                                entered_exe=bedtobigbed_exe)
        else:
            # Cloud mode; check S3 for genome if necessary
            assert s3_ansible is not None
            if assembly == 'hg19'
                base.index_archive = 's3://rail-emr/index/hg19_UCSC.tar.gz'
            else:
                if not Url(assembly).is_s3:
                    base.errors.append('Bowtie index archive must be on S3 '
                                       ' in cloud ("--cloud") mode.')
                elif not s3_ansible.exists(assembly):
                    base.errors.append('Bowtie index archive was not found '
                                       'on S3 at "{0}".'.format(assembly))
        # Assume bowtie2 args are kosher for now
        base.bowtie2_args = bowtie2_args
        if not (isinstance(genome_partition_length, int) and
                genome_partition_length > 0):
            base.errors.append('Genome partition length '
                               '(--genome-partition-length) must be an '
                               'integer > 0, but {0} was entered.'.format(
                                                        genome_partition_length
                                                    ))
        base.genome_partition_length = genome_partition_length
        if not (isinstance(min_readlet_size, int) and min_readlet_size > 0):
            base.errors.append('Minimum readlet size (--min-readlet-size) '
                               'must be an integer > 0, but '
                               '{0} was entered.'.format(min_readlet_size))
        base.min_readlet_size = min_readlet_size
        if not (isinstance(max_readlet_size, int) and max_readlet_size
                >= min_readlet_size):
            base.errors.append('Maximum readlet size (--max-readlet-size) '
                               'must be an integer >= minimum readlet size '
                               '(--min-readlet-size) = '
                               '{0}, but {1} was entered.'.format(
                                                    base.min_readlet_size,
                                                    max_readlet_size
                                                ))
        base.max_readlet_size = max_readlet_size
        if not (isinstance(readlet_interval, int) and readlet_interval
                > 0):
            base.errors.append('Readlet interval (--readlet-interval) '
                               'must be an integer > 0, '
                               'but {0} was entered.'.format(
                                                    readlet_interval
                                                ))
        base.readlet_interval = readlet_interval
        if not (isinstance(cap_size_multiplier, float) and cap_size_multiplier
                > 1):
            base.errors.append('Cap size multiplier (--cap-size-multiplier) '
                               'must be > 1, '
                               'but {0} was entered.'.format(
                                                    cap_size_multiplier
                                                ))
        base.cap_size_multiplier = cap_size_multiplier
        if not (isinstance(min_intron_size, int) and min_intron_size > 0):
            base.errors.append('Minimum intron size (--min-intron-size) '
                               'must be an integer > 0, but '
                               '{0} was entered.'.format(min_intron_size))
        base.min_intron_size = min_intron_size
        if not (isinstance(max_intron_size, int) and max_intron_size
                >= min_intron_size):
            base.errors.append('Maximum intron size (--max-intron-size) '
                               'must be an integer >= minimum intron size '
                               '(--min-readlet-size) = '
                               '{0}, but {1} was entered.'.format(
                                                    base.min_intron_size,
                                                    max_intron_size
                                                ))
        base.max_intron_size = max_intron_size
        if not (isinstance(min_exon_size, int) and min_exon_size > 0):
            base.errors.append('Minimum exon size (--min-exon-size) '
                               'must be an integer > 0, but '
                               '{0} was entered.'.format(min_exon_size))
        base.min_exon_size = min_exon_size
        if not (isinstance(motif_search_window_size, int) and 
                    motif_search_window_size >= 0):
            base.errors.append('Motif search window size '
                               '(--motif-search-window-size) must be an '
                               'integer > 0, but {0} was entered.'.format(
                                                    motif_search_window_size
                                                ))
        if not (isinstance(motif_radius, int) and
                    motif_radius >= 0):
            base.errors.append('Motif radius (--motif-radius) must be an '
                               'integer >= 0, but {0} was entered.'.format(
                                                    motif_radius
                                                ))
        base.motif_radius = motif_radius
        if not (isinstance(normalize_percentile, float) and
                    0 <= normalize_percentile <= 1):
            base.errors.append('Normalization percentile '
                               '(--normalize-percentile) must on the '
                               'interval [0, 1], but {0} was entered'.format(
                                                    normalize_percentile
                                                ))
        base.normalize_percentile = normalize_percentile
        base.do_not_output_bam_by_chr = do_not_output_bam_by_chr
        base.output_sam = output_sam
        base.bam_basename = bam_basename
        base.bed_basename = bed_basename

    @staticmethod
    def add_args(parser, usage=None, local=False):
        """ usage: argparse.SUPPRESS if advanced options should be suppressed;
                else None
        """
        if local:
            parser.add_argument(
                '--bowtie1-exe', type=str, required=False,
                default=None,
                help=('Path to Bowtie 1 executable. This can be left out if '
                      '"bowtie" is in PATH and is executable.')
            )
            parser.add_argument(
                '-1', '--bowtie1-idx', type=str, required=True,
                help='Path to Bowtie 1 index. Include basename.'
            )
            parser.add_argument(
                '--bowtie2-exe', type=str, required=False,
                default=None,
                help=('Path to Bowtie 2 executable. This can be left out if '
                      '"bowtie2" is in PATH and is executable.')
            )
            parser.add_argument(
                '-2', '--bowtie2-idx', type=str, required=True,
                help='Path to Bowtie 2 index. Include basename.'
            )
            parser.add_argument(
                '--bowtie2_args', type=str, required=False,
                default='',
                help=('Additional arguments to pass to Bowtie 2, which is '
                      'used to obtain final end-to-end alignments and spliced '
                      'alignments in output SAM/BAM files.')
            )
            parser.add_argument(
                '--samtools-exe', type=str, required=False,
                default=None,
                help=('Path to SAMTools executable. This can be left out if '
                      '"samtools" is in PATH and is executable.')
            )
            parser.add_argument(
                '--bedtobigbed-exe', type=str, required=False,
                default=None,
                help=('Path to BedToBigBed executable. This can be left out '
                      'if "bedToBigBed" is in PATH and is executable.')
            )
        else:
            parser.add_argument(
                '--assembly', type=str, required=False,
                default='hg19',
                help=('One of the following:\n'
                      '   a) An assembly identifier for a supported '
                      'organism. Currently supported: {"hg19"}.\n'
                      '   b) Accessible S3 path to a Rail archive. A Rail '
                      'archive is a tar.gz composed a folder titled "index", '
                      'that contains BOTH Bowtie 1 and Bowtie 2 index files '
                      '(12 in total), all with the basename "genome".')
            )
        parser.add_argument(
            '--genome-partition-length', type=int, required=False,
            usage=usage,
            default=5000,
            help=('Smallest unit of a genome (in nucleotides) addressable by '
                  'a single task when computing coverage from exon '
                  'differentials. Making this parameter too small (~hundreds '
                  'of nucleotides) can bloat intermediate files, but making '
                  'it too large (~the size of a chromosome) could compromise '
                  'the scalability of the pipeline.')
        )
        parser.add_argument(
            '--max-readlet-size', type=int, required=False,
            default=25,
            usage=usage,
            help=('Maximum size of a given segment from a read that is 1) '
                  'mapped to the genome using Bowtie 1 when searching for '
                  'introns; and 2) mapped to a set of transcriptome elements '
                  'using Bowtie 2 when finalizing spliced alignments.'
                  'Decreasing the value of this parameter may increase recall '
                  'of introns while compromising precision. For human-size '
                  'genomes, values between 20 and 25 are recommended.')
        )
        parser.add_argument(
            '--min-readlet-size', type=int, required=False,
            default=15,
            usage=usage,
            help=('Minimum size of a given "capping readlet" -- that is, a '
                  'read segment whose end coincides with a read end -- that '
                  'is 1) mapped to the genome using Bowtie 1 when searching '
                  'for introns; and 2) mapped to a set of transcriptome '
                  'elements using Bowtie 2 when finalizing spliced '
                  'alignments. Decreasing the value of this parameter may '
                  'increase recall of rare introns overlapped towards the '
                  'ends of a few reads in a sample while compromising '
                  'precision.')
        )
        parser.add_argument(
            '--readlet-interval', type=int, required=False,
            default=4,
            usage=usage,
            help=('Distance (in nucleotides) between overlapping readlets '
                  'mapped to genome and set of transcriptome elements. '
                  'Decreasing this parameter may increase sensitivity while '
                  'increasing the time the pipeline takes.')
        )
        parser.add_argument(
            '--cap-size-multiplier', type=float, required=False,
            default=1.2,
            usage=usage,
            help=('Successive capping readlets on a given end of a read are '
                  'increased in size exponentially with this base.')
        )
        parser.add_argument(
            '--max-intron-size', type=int, required=False,
            default=500000,
            help=('Introns spanning more than this many nucleotides are '
                  'automatically filtered out.')
        )
        parser.add_argument(
            '--min-intron-size', type=int, required=False,
            default=10,
            help=('Introns spanning fewer than this many nucleotides are '
                  'automatically filtered out.')
        )
        parser.add_argument(
            '--min-exon-size', type=int, required=False,
            default=9,
            help=('The aligner will not be sensitive to exons smaller than '
                  'this size.')
        )
        parser.add_argument(
            '--motif-search-window-size', type=int, required=False,
            default=1000,
            usage=usage,
            help=('Size of window in which to search for exons of size '
                  '"--min-exon-size" capped by appropriate donor/acceptor '
                  'motifs when inferring intron positions.')
        )
        parser.add_argument(
            '--motif-radius', type=int, required=False,
            default=5,
            usage=usage,
            help=('Number of nucleotides of ostensible intron ends within '
                  'which to search for a donor/acceptor motif.')
        )
        parser.add_argument(
            '--normalize-percentile', type=float, required=False,
            default=0.75,
            help=('Percentile used for computing normalization factors for '
                  'sample coverage.')
        )
        parser.add_argument(
            '--do-not-output-bam-by-chr', action='store_const', const=True,
            default=False,
            help=('Places alignments for all chromosomes in a single file '
                  'rather than dividing them up by reference name.')
        )
        parser.add_argument(
            '--output-sam', action='store_const', const=True,
            default=False,
            help='Outputs SAM files rather than BAM files.'
        )
        parser.add_argument(
            '--bam-basename', type=str, required=False,
            default='alignments',
            help='Basename to use for BAM output files.'
        )
        parser.add_argument(
            '--bed-basename', type=str, required=False,
            default='',
            help=('Basename to use for BED output files; there is an output '
                  'for each of insertions, deletions, and introns.')
        )

    @staticmethod
    def pseudosteps(input_dir, intermediate_dir, output_dir, step_dir):


class RailRnaLocalPreprocessJson:
    """ Constructs JSON for local mode + preprocess job flow. """
    def __init__(self, manifest, output_dir, intermediate_dir='./intermediate',
        force=False, aws_exe=None, profile='default', region='us-east-1',
        verbose=False, nucleotides_per_input=8000000, gzip_input=True,
        num_processes=1, keep_intermediates=False):
        base = RailRnaErrors(manifest, output_dir, 
            intermediate_dir=intermediate_dir,
            force=force, aws_exe=aws_exe, profile=profile,
            region=region, verbose=verbose)
        RailRnaLocal(base, check_manifest=True, num_processes=num_processes,
            keep_intermediates=keep_intermediates)
        RailRnaPreprocess(base, nucleotides_per_input=nucleotides_per_input,
            gzip_input=gzip_input)
        if base.errors:
            raise RuntimeError(
                    '\n'.join(
                            ['%d) %s' % (i, error) for i, error in errors]
                        )
                )
        self._json_serial = {}
        self._json_serial['Steps'] = []
        steps = self._json_serial['Steps']
        step_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 
                                                    'steps'))
        steps.append(
            step(name='Preprocess input reads',
                    tasks=0,
                    inputs=base.manifest,
                    output=os.path.join(base_intermediate_dir, 'preprocess'),
                    mapper=' '.join([sys.executable, 
                                        os.path.join(step_dir, 
                                            'preprocess.py'),
                                        '--nucs-per-file=%d'
                                        % base.nucleotides_per_input,
                                        '--gzip-output'
                                        if base.gzip_input else '',
                                        '--ignore-first-token',
                                        '--push=',
                                        base.output_dir])
                )
            )
    
    @property
    def json_serial(self):
        return self._json_serial
