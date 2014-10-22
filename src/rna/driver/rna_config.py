#!/usr/bin/env python
"""
rna_config.py
Part of Rail-RNA

Contains classes that perform error checking and generate JSON configuration
output for Rail-RNA. These configurations are parsable by Dooplicity's
emr_simulator.py and emr_runner.py.

Class structure is designed so only those arguments relevant to modes/job flows
are included.

The descriptions of command-line arguments contained here assume that a
calling script has the command-line options "prep", "align", "go", "local",
and "elastic".
"""

import os
base_path = os.path.abspath(
                    os.path.dirname(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__)))
                    )
                )
utils_path = os.path.join(base_path, 'rna', 'utils')
import site
site.addsitedir(utils_path)
site.addsitedir(base_path)
import dooplicity.ansibles as ab
import tempfile
import shutil
from dooplicity.tools import path_join, is_exe, which
from version import version_number
import sys
from argparse import SUPPRESS
import subprocess
import time
from traceback import format_exc
from collections import defaultdict

'''These are placed here for convenience; their locations may change
on EMR depending on bootstraps.'''
_hadoop_streaming_jar = '/home/hadoop/contrib/streaming/hadoop-streaming.jar'
_multiple_files_jar = '/mnt/lib/multiple-files.jar'
_relevant_elephant_jar = '/mnt/lib/relevant-elephant.jar'
_hadoop_lzo_jar = ('/home/hadoop/.versions/2.4.0/share/hadoop'
                   '/common/lib/hadoop-lzo.jar')
_s3distcp_jar = '/home/hadoop/lib/emr-s3distcp-1.0.jar'
_hdfs_temp_dir = 'hdfs:///railtemp'
_base_combine_split_size = 268435456 # 250 MB
_elastic_bowtie1_idx = '/mnt/index/genome'
_elastic_bowtie2_idx = '/mnt/index/genome'
_elastic_bedgraphtobigwig_exe ='/mnt/bin/bedGraphToBigWig'
_elastic_samtools_exe = 'samtools'
_elastic_bowtie1_exe = 'bowtie'
_elastic_bowtie2_exe = 'bowtie2'
_elastic_bowtie1_build_exe = 'bowtie-build'
_elastic_bowtie2_build_exe = 'bowtie2-build'

# Decide Python executable
if 'pypy 2.' in sys.version.lower():
    # Executable has the user's desired version of PyPy
    _warning_message = 'Launching Dooplicity runner with PyPy...'
    _executable = sys.executable
else:
    _pypy_exe = which('pypy')
    _print_warning = False
    if _pypy_exe is not None:
        try:
            if 'pypy 2.' in \
                subprocess.check_output(
                        [_pypy_exe, '--version'],
                        stderr=subprocess.STDOUT
                    ).lower():
                _executable = _pypy_exe
            else:
                _executable = sys.executable
                _print_warning = True
        except Exception as e:
            _executable = sys.executable
            _print_warning = True
    else:
        _print_warning = True
    if _print_warning:
        _warning_message = ('WARNING: PyPy 2.x not found. '
            'Installation is recommended to optimize performance '
            'of Rail-RNA in "local" or "parallel" mode. If it is installed, '
            'make sure the "pypy" executable is in PATH, or use '
            'it to execute Rail-RNA via "[path to \'pypy\'] '
            '[path to \'rail-rna\']".')
        _executable = sys.executable
    else:
        _warning_message = 'Launching Dooplicity runner with PyPy...'

def ready_engines(rc, base):
    """ Prepares engines for checks. 

        rc: IPython Client object
        base: instance of RailRnaErrors

        No return value.
    """
    rc[:].execute('import socket').get()
    rc[:].execute('import subprocess').get()
    engine_to_hostnames = rc[:].apply(socket.get_hostnames).get_dict()
    hostname_to_engines = defaultdict(set)
    for engine in hostnames:
        hostname_to_engines[engine_to_hostnames[engine]].add(engine)
    # Distribute Rail-RNA to nodes
    
    while not asyncresult.ready():
        time.sleep(1e-1)
    if not asyncresult.successful():
        asyncresult.get()
    asyncresult = rc[:]
    asyncresult = rc[:].execute('import os')
    while not asyncresult.ready():
        time.sleep(1e-1)
    if not asyncresult.successful():
        asyncresult.get()

def step(name, inputs, output,
    mapper='org.apache.hadoop.mapred.lib.IdentityMapper',
    reducer='org.apache.hadoop.mapred.lib.IdentityReducer', 
    action_on_failure='TERMINATE_JOB_FLOW',
    jar=_hadoop_streaming_jar,
    tasks=0, partitioner_options=None, key_fields=None, archives=None,
    multiple_outputs=False, inputformat=None, extra_args=[]):
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
        extra_args: extra '-D' args

        Return value: step dictionary
    """
    to_return = {
        'Name' : name,
        'ActionOnFailure' : action_on_failure,
        'HadoopJarStep' : {
            'Jar' : jar,
            'Args' : []
        }
    }
    to_return['HadoopJarStep']['Args'].extend(
            ['-D', 'mapreduce.job.reduces=%d' % tasks]
        )
    if partitioner_options is not None and key_fields is not None:
        to_return['HadoopJarStep']['Args'].extend([
                '-D', 'mapreduce.partition.keypartitioner.options=-%s'
                            % partitioner_options,
                '-D', 'stream.num.map.output.key.fields=%d' % key_fields
            ])
    for extra_arg in extra_args:
        to_return['HadoopJarStep']['Args'].extend(
            ['-D', extra_arg]
        )
    # Add libjar for splittable LZO
    to_return['HadoopJarStep']['Args'].extend(
            ['-libjars', _relevant_elephant_jar]
        )
    if multiple_outputs:
        to_return['HadoopJarStep']['Args'][-1] \
            +=  (',%s' % _multiple_files_jar)
    if archives is not None:
        to_return['HadoopJarStep']['Args'].extend([
                '-archives', archives
            ])
    to_return['HadoopJarStep']['Args'].extend([
            '-partitioner',
            'org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner',
        ])
    to_return['HadoopJarStep']['Args'].extend([
                '-input', ','.join([an_input.strip() for an_input in inputs])
            ])
    to_return['HadoopJarStep']['Args'].extend([
            '-output', output,
            '-mapper', mapper,
            '-reducer', reducer
        ])
    if multiple_outputs:
        to_return['HadoopJarStep']['Args'].extend([
                '-outputformat', 'edu.jhu.cs.MultipleOutputFormat'
            ])
    if inputformat is not None:
        to_return['HadoopJarStep']['Args'].extend([
                '-inputformat', inputformat
            ])
    else:
        '''Always use splittable LZO; it's deprecated because hadoop-streaming
        uses the old mapred API.'''
        to_return['HadoopJarStep']['Args'].extend([
                '-inputformat',
                'com.twitter.elephantbird.mapred.input'
                '.DeprecatedCombineLzoTextInputFormat'
            ])
    return to_return

# TODO: Flesh out specification of protostep and migrate to Dooplicity
def steps(protosteps, action_on_failure, jar, step_dir, 
            reducer_count, intermediate_dir, extra_args=[], unix=False,
            no_consistent_view=False):
    """ Turns list with "protosteps" into well-formed StepConfig list.

        A protostep looks like this:

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
                'taskx' : number of tasks per reducer or None if total number
                    of tasks should be 1
                'inputformat' : input format; present only if necessary
                'archives' : archives parameter; present only if necessary
                'multiple_outputs' : key that's present iff there are multiple
                    outputs
                'index_output' : key that's present iff output LZOs should be
                    indexed after step; applicable only in Hadoop modes
                'direct_copy' : if output directory is s3, copy outputs there 
                    directly; do not use hdfs
                'extra_args' : list of '-D' args
            }

        protosteps: array of protosteps
        action_on_failure: action on failure to take
        jar: path to Hadoop Streaming jar
        step_dir: where to find Python scripts for steps
        reducer_count: number of reducers; determines number of tasks
        unix: performs UNIX-like path joins; also inserts pypy in for
            executable since unix=True only on EMR
        no_consistent_view: True iff consistent view should be switched off;
            adds s3distcp commands when there are multiple outputs

        Return value: list of StepConfigs (see Elastic MapReduce API docs)
    """
    '''CombineFileInputFormat outputs file offsets as keys; kill them here
    in an "identity" mapper. This could also be done by overriding an
    appropriate method from class in Java, but the overhead incurred doing
    it this way should be small.'''
    true_steps = []
    for protostep in protosteps:
        assert ('keys' in protostep and 'part' in protostep) or \
                ('keys' not in protostep and 'part' not in protostep)
        # Can't copy directly to S3 if there are multiple outputs ... yet!
        assert not (('direct_copy' in protostep) and ('multiple_outputs'
                        in protostep))
        identity_mapper = ('cut -f 2-' if unix else 'cat')
        final_output = (path_join(unix, intermediate_dir,
                                        protostep['output'])
                        if 'no_output_prefix' not in
                        protostep else protostep['output'])
        final_output_url = ab.Url(final_output)
        if (not ('direct_copy' in protostep) and unix
            and final_output_url.is_s3 and no_consistent_view):
            intermediate_output = _hdfs_temp_dir + final_output_url.suffix[1:]
        else:
            intermediate_output = final_output
        true_steps.append(step(
                name=protostep['name'],
                inputs=([path_join(unix, intermediate_dir,
                            an_input) for an_input in
                            protostep['inputs']]
                        if 'no_input_prefix' not in
                        protostep else protostep['inputs']),
                output=intermediate_output,
                mapper=' '.join(['pypy' if unix
                        else _executable, 
                        path_join(unix, step_dir,
                                        protostep['run'])])
                        if 'keys' not in protostep else identity_mapper,
                reducer=' '.join(['pypy' if unix
                        else _executable, 
                        path_join(unix, step_dir,
                                        protostep['run'])]) 
                        if 'keys' in protostep else 'cat',
                action_on_failure=action_on_failure,
                jar=jar,
                tasks=((reducer_count * protostep['taskx'] * 10 / 8)
                        if protostep['taskx'] is not None
                        else 1),
                partitioner_options=(protostep['part']
                    if 'part' in protostep else None),
                key_fields=(protostep['keys']
                    if 'keys' in protostep else None),
                archives=(protostep['archives']
                    if 'archives' in protostep else None),
                multiple_outputs=(True if 'multiple_outputs'
                        in protostep else False
                    ),
                inputformat=(protostep['inputformat']
                    if 'inputformat' in protostep else None),
                extra_args=(protostep['extra_args']
                    if 'extra_args' in protostep else [])
            )
        )
        if unix and 'index_output' in protostep:
            # Index inputs before subsequent step begins
            true_steps.append(
                    {
                        'Name' : ('Index output of "'
                                    + protostep['name'] + '"'),
                        'ActionOnFailure' : action_on_failure,
                        'HadoopJarStep' : {
                            'Jar' : _hadoop_lzo_jar,
                            'Args' : ['com.hadoop.compression.lzo'
                                      '.DistributedLzoIndexer',
                                      intermediate_output]
                        }
                    }
                )
        if (not ('direct_copy' in protostep) and unix
            and final_output_url.is_s3 and no_consistent_view):
            # s3distcp intermediates over
            true_steps.append(
                    {
                        'Name' : ('Copy output of "'
                                    + protostep['name'] + '" to S3'),
                        'ActionOnFailure' : action_on_failure,
                        'HadoopJarStep' : {
                            'Jar' : _s3distcp_jar,
                            'Args' : ['--src', intermediate_output,
                                      '--dest', final_output,
                                      '--deleteOnSuccess']
                        }
                    }
                )
    return true_steps

class RailRnaErrors(object):
    """ Holds accumulated errors in Rail-RNA's input parameters.

        Checks only those parameters common to all modes/job flows.
    """
    def __init__(self, manifest, output_dir,
            intermediate_dir='./intermediate', force=False, aws_exe=None,
            profile='default', region='us-east-1', verbose=False,
            curl_exe=None
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
        self.profile = profile

    def check_s3(self, reason=None, is_exe=None, which=None):
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
        if not is_exe:
            is_exe = globals()['is_exe']
        if not which:
            which = globals()['which']
        original_errors_size = len(self.errors)
        if self.aws_exe is None:
            self.aws_exe = 'aws'
            if not which(self.aws_exe):
                self.errors.append(('The AWS CLI executable '
                                    'was not found. Make sure that the '
                                    'executable is in PATH, or specify the '
                                    'location of the executable with '
                                    '--aws.'))
        elif not is_exe(self.aws_exe):
            self.errors.append(('The AWS CLI executable (--aws) '
                                '"{0}" is either not present or not '
                                'executable.').format(aws_exe))
        self._aws_access_key_id = None
        self._aws_secret_access_key = None
        if self.profile == 'default':
            # Search environment variables for keys first if profile is default
            try:
                self._aws_access_key_id = os.environ['AWS_ACCESS_KEY_ID']
                self._aws_secret_access_key \
                    = os.environ['AWS_SECRET_ACCESS_KEY']
            except KeyError:
                to_search = '[default]'
            else:
                to_search = None
            try:
                # Also grab region
                self.region = os.environ['AWS_DEFAULT_REGION']
            except KeyError:
                pass
        else:
            to_search = '[profile ' + self.profile + ']'
        # Now search AWS CLI config file for the right profile
        if to_search is not None:
            cred_file = os.path.join(os.environ['HOME'], '.aws', 'cred')
            if os.path.exists(cred_file):
                # "credentials" file takes precedence over "config" file
                config_file = cred_file
            else:
                config_file = os.path.join(os.environ['HOME'], '.aws',
                                            'config')
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
                                    'the profile (--profile) is set to '
                                    '"default" (its default value).\n\n'
                                    'b) The file ".aws/config" or '
                                    '".aws/credentials" exists in your '
                                    'home directory with a valid profile. '
                                    'To set this file up, run "aws --config" '
                                    'after installing the AWS CLI.')
                                )
        if len(self.errors) != original_errors_size:
            if reason:
                raise RuntimeError((('\n'.join(['%d) %s' % (i+1, error)
                                    for i, error
                                    in enumerate(self.errors)]) 
                                    if len(self.errors) > 1
                                    else self.errors[0]) + 
                                    '\n\nNote that the AWS CLI is needed '
                                    'because {0}. If all dependence on S3 in '
                                    'the pipeline is removed, the AWS CLI '
                                    'need not be installed.').format(reason))
            else:
                raise RuntimeError((('\n'.join(['%d) %s' % (i+1, error)
                                    for i, error
                                    in enumerate(self.errors)])
                                    if len(self.errors) > 1
                                    else self.errors[0]) + 
                                    '\n\nIf all dependence on S3 in the '
                                    'pipeline is removed, the AWS CLI need '
                                    'not be installed.'))
        self.checked_programs.add('AWS CLI')

    def check_program(self, exe, program_name, parameter,
                        entered_exe=None, reason=None,
                        is_exe=None, which=None):
        """ Checks if program in PATH or if user specified it properly.

            Errors are added to self.errors.

            exe: executable to search for
            program name: name of program
            parameter: corresponding command line parameter
                (e.g., --bowtie)
            entered_exe: None if the user didn't enter an executable; otherwise
                whatever the user entered
            reason: FOR CURL ONLY: raise RuntimeError _immediately_ if Curl
                not found but needed

            No return value.
        """
        if not is_exe:
            is_exe = globals()['is_exe']
        if not which:
            which = globals()['which']
        original_errors_size = len(self.errors)
        if entered_exe is None:
            if not which(exe):
                self.errors.append(
                        ('The executable "{0}" for {1} was either not found '
                         'in PATH or is not executable. Check that the '
                         'program is installed properly and executable; then '
                         'either add the executable to PATH or specify it '
                         'directly with {2}.').format(exe, program_name,
                                                            parameter)
                    )
            to_return = exe
        elif not is_exe(entered_exe):
            self.errors.append(
                    ('The executable "{0}" entered for {1} via {2} was '
                     'either not found or is not executable.').format(exe,
                                                                program_name,
                                                                parameter)
                )
            to_return = entered_exe
        else:
            to_return = entered_exe
        if original_errors_size != len(self.errors) and reason:
            raise RuntimeError((('\n'.join(['%d) %s' % (i+1, error)
                                for i, error
                                in enumerate(self.errors)])
                                if len(self.errors) > 1 else self.errors[0]) + 
                                '\n\nNote that Curl is needed because {0}.'
                                ' If all dependence on web resources is '
                                'removed from the pipeline, Curl need '
                                'not be installed.').format(reason))
        self.checked_programs.add(program_name)
        return to_return

    @staticmethod
    def add_args(general_parser, exec_parser, required_parser):
        exec_parser.add_argument(
            '--aws', type=str, required=False, metavar='<exe>',
            default=None,
            help='path to AWS CLI executable (def: aws)'
        )
        exec_parser.add_argument(
            '--curl', type=str, required=False, metavar='<exe>',
            default=None,
            help='path to Curl executable (def: curl)'
        )
        general_parser.add_argument(
            '--profile', type=str, required=False, metavar='<str>',
            default='default',
            help='AWS CLI profile (def: env vars, then "default")'
        )
        general_parser.add_argument(
            '-f', '--force', action='store_const', const=True,
            default=False,
            help='overwrite output directory if it exists'
        )
        general_parser.add_argument(
            '--verbose', action='store_const', const=True,
            default=False,
            help='write extra debugging statements to stderr'
        )
        required_parser.add_argument(
            '-m', '--manifest', type=str, required=True, metavar='<file>',
            help='Myrna-style manifest file; Google "Myrna manifest" for ' \
                 'help'
        )
        '''--region's help looks different from mode to mode; don't include it
        here.'''

def ipython_client(ipython_profile=None, ipcontroller_json=None):
    """ Performs checks on IPython engines and returns IPython Client object.

        Also checks that IPython is installed/configured properly and raises
        exception _immediately_ if it's not; then prints
        engine detection message.

        ipython_profile: IPython parallel profile, if specified; otherwise None
        ipcontroller_json: IP Controller json file, if specified; else None

        No return value.
    """
    errors = []
    try:
        from IPython.parallel import Client
    except ImportError:
        errors.append(
                   'IPython is required to run Rail-RNA in '
                   '"parallel" mode. Visit ipython.org to '
                   'download it, or simply download the Anaconda '
                   'distribution of Python at '
                   'https://store.continuum.io/cshop/anaconda/; it\'s '
                   'easy to install and comes with IPython and '
                   'several other useful packages.'
                )
    if ipython_profile:
        try:
            rc = Client(profile=ipy_profile)
        except ValueError:
            errors.append(
                    'Cluster configuration profile "%s" was not '
                    'found.' % ipy_profile
                )
    elif ipcontroller_json:
        try:
            rc = Client(ipcontroller_json)
        except IOError:
            errors.append(
                    'Cannot find connection information JSON file %s.'
                    % ipcontroller_json
                )
    else:
        try:
            rc = Client()
        except IOError:
            errors.append(
                    'Cannot find ipcontroller-client.json. Ensure '
                    'that IPython controller and engines are running.'
                    ' If controller is running on a remote machine, '
                    'copy the ipcontroller-client.json file from there '
                    'to a local directory; then rerun this script '
                    'specifying the local path to '
                    'ipcontroller-client.json with the '
                    '--ipcontroller-json command-line parameter.'
                )
    if errors:
        raise RuntimeError(
                    '\n'.join(
                            ['%d) %s' % (i+1, error) for i, error
                                in enumerate(errors)]
                        ) if len(errors) > 1 else errors[0]
                )
    if not rc.ids:
        raise RuntimeError(
                'An IPython controller is running, but no engines are '
                'connected to it. Engines must be connected to an IPython '
                'controller when running Rail-RNA in "parallel" mode.'
            )
    engine_detect_message = ('Detected %d running IPython engines.' 
                                % len(rc.ids))
    print >>sys.stderr, engine_detect_message
    if not sys.stderr.isatty():
        # So the user sees it too
        print engine_detect_message
    return rc

class RailRnaLocal(object):
    """ Checks local- or parallel-mode JSON for programs and input parameters.

        Subsumes only those parameters relevant to local mode. Adds errors
        to base instance of RailRnaErrors.
    """
    def __init__(self, base, check_manifest=False,
                    num_processes=1, keep_intermediates=False,
                    gzip_intermediates=False, gzip_level=3, parallel=False,
                    local=True, scratch=None, ansible=None):
        """ base: instance of RailRnaErrors """
        # Initialize ansible for easy checks
        if not ansible:
            ansible = ab.Ansible()
        if not ab.Url(base.intermediate_dir).is_local:
            base.errors.append(('Intermediate directory must be in locally '
                                'accessible filesystem when running Rail-RNA '
                                'in "local" or "parallel" mode, but {0} was '
                                'entered.').format(
                                        base.intermediate_dir
                                    ))
        else:
            base.intermediate_dir = os.path.abspath(base.intermediate_dir)
        output_dir_url = ab.Url(base.output_dir)
        if output_dir_url.is_curlable:
            base.errors.append(('Output directory must be in locally '
                                'accessible filesystem or on S3 '
                                'when running Rail-RNA in "local" '
                                'or "parallel" mode, '
                                'but {0} was entered.').format(
                                        base.output_dir
                                    ))
        elif output_dir_url.is_s3 and 'AWS CLI' not in base.checked_programs:
            base.check_s3(reason='the output directory is on S3',
                            is_exe=is_exe,
                            which=which)
            # Change ansible params
            ansible.aws_exe = base.aws_exe
            ansible.profile = base.profile
        base.check_s3_on_engines = None
        base.check_curl_on_engines = None
        if not parallel:
            if output_dir_url.is_local \
                and os.path.exists(output_dir_url.to_url()):
                if not base.force:
                    base.errors.append(('Output directory {0} exists, '
                                        'and --force was not invoked to '
                                        'permit overwriting it.').format(
                                                base.output_dir)
                                            )
                else:
                    try:
                        shutil.rmtree(base.output_dir)
                    except OSError:
                        try:
                            os.remove(base.output_dir)
                        except OSError:
                            pass
                    base.output_dir = os.path.abspath(base.output_dir)
            elif output_dir_url.is_s3 \
                and ansible.s3_ansible.is_dir(base.output_dir):
                if not base.force:
                    base.errors.append(('Output directory {0} exists on S3, '
                                        'and --force was not invoked to '
                                        'permit overwriting it.').format(
                                                base_output_dir)
                                            )
                else:
                    ansible.s3ansible.remove_dir(base.output_dir)
            # Check manifest; download it if necessary
            manifest_url = ab.Url(base.manifest)
            if manifest_url.is_s3 and 'AWS CLI' not in base.checked_programs:
                base.check_s3(reason='the manifest file is on S3',
                                is_exe=is_exe,
                                which=which)
                # Change ansible params
                ansible.aws_exe = base.aws_exe
                ansible.profile = base.profile
            elif manifest_url.is_curlable \
                and 'Curl' not in base.checked_programs:
                base.curl_exe = base.check_program('curl', 'Curl', '--curl',
                                    entered_exe=base.curl_exe,
                                    reason='the manifest file is on the web',
                                    is_exe=is_exe,
                                    which=which)
                ansible.curl_exe = base.curl_exe
            if not ansible.exists(manifest_url.to_url()):
                base.errors.append(('Manifest file (--manifest) {0} '
                                    'does not exist. Check the URL and '
                                    'try again.').format(base.manifest))
            else:
                if not manifest_url.is_local:
                    '''Download/check manifest only if not an IPython engine
                    (not parallel).'''
                    base.manifest_dir = base.intermediate_dir
                    import atexit
                    from tempdel import remove_temporary_directories
                    atexit.register(remove_temporary_directories,
                                        [base.manifest_dir])
                    base.manifest = os.path.join(base.manifest_dir, 'MANIFEST')
                    ansible.get(manifest_url, destination=base.manifest)
                base.manifest = os.path.abspath(base.manifest)
                files_to_check = []
                base.sample_count = 0
                with open(base.manifest) as manifest_stream:
                    for line in manifest_stream:
                        if line[0] == '#' or not line.strip(): continue
                        base.sample_count += 1
                        tokens = line.strip().split('\t')
                        check_sample_label = True
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
                                                        manifest_url.to_url(),
                                                        line
                                                    ))
                            check_sample_label = False
                        if check_sample_label and tokens[-1].count('-') != 2:
                            line = line.strip()
                            base.errors.append(('The following line from the '
                                                'manifest file {0} '
                                                'has an invalid sample label: '
                                                '\n{1}\nA valid sample label '
                                                'takes the following form:\n'
                                                '<Group ID>-<BioRep ID>-'
                                                '<TechRep ID>'
                                                ).format(
                                                        manifest_url.to_url(),
                                                        line
                                                    ))
                if files_to_check:
                    if check_manifest:
                        # Check files in manifest only if in preprocess flow
                        file_count = len(files_to_check)
                        for k, filename in enumerate(files_to_check):
                            if sys.stdout.isatty():
                                sys.stdout.write(
                                        '\r\x1b[KChecking that file %d/%d '
                                        'from manifest file exists...' % (
                                                                    k+1,
                                                                    file_count
                                                                )
                                    )
                                sys.stdout.flush()
                            filename_url = ab.Url(filename)
                            if filename_url.is_s3 \
                                and 'AWS CLI' not in base.checked_programs:
                                    if local:
                                        base.check_s3(reason=(
                                                      'at least one sample '
                                                      'FASTA/FASTQ from the '
                                                      'manifest file is on '
                                                      'S3'),
                                                      is_exe=is_exe,
                                                      which=which
                                                    )
                                    base.check_s3_on_engines = (
                                                      'at least one sample '
                                                      'FASTA/FASTQ from the '
                                                      'manifest file is on '
                                                      'S3'
                                                    )
                                    # Change ansible params
                                    ansible.aws_exe = base.aws_exe
                                    ansible.profile = base.profile
                            elif filename_url.is_curlable \
                                and 'Curl' not in base.checked_programs:
                                if local:
                                    base.curl_exe = base.check_program('curl',
                                                    'Curl',
                                                    '--curl',
                                                    entered_exe=base.curl_exe,
                                                    reason=(
                                                      'at least one sample '
                                                      'FASTA/FASTQ from the '
                                                      'manifest file is on '
                                                      'the web'),
                                                    is_exe=is_exe,
                                                    which=which
                                                )
                                base.check_curl_on_engines = (
                                                      'at least one sample '
                                                      'FASTA/FASTQ from the '
                                                      'manifest file is on '
                                                      'the web'
                                                    )
                                ansible.curl_exe = base.curl_exe
                            if not ansible.exists(filename):
                                base.errors.append((
                                                    'The file {0} from the '
                                                    'manifest file {1} does '
                                                    'not exist. Check the URL '
                                                    'and try again.').format(
                                                        filename,
                                                        manifest_url.to_url()
                                                ))
                        if sys.stdout.isatty():
                            sys.stdout.write(
                                    '\r\x1b[KChecked all files listed in '
                                    'manifest file.\n'
                                )
                            sys.stdout.flush()
                else:
                    base.errors.append(('Manifest file (--manifest) {0} '
                                        'has no valid lines.').format(
                                                        manifest_url.to_url()
                                                    ))
            from multiprocessing import cpu_count
            if num_processes:
                if not (isinstance(num_processes, int)
                                        and num_processes >= 1):
                    base.errors.append('Number of processes (--num-processes) '
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
                    so Facebook tab in user's browser won't go all
                    unresponsive.'''
                    base.num_processes -= 1
            if gzip_intermediates:
                if not (isinstance(gzip_level, int)
                                        and 9 >= gzip_level >= 1):
                    base.errors.append('Gzip level (--gzip-level) '
                                       'must be an integer between 1 and 9, '
                                       'but {0} was entered.'.format(
                                                        gzip_level
                                                    ))
        if scratch:
            if not os.path.exists(scratch):
                try:
                    os.makedirs(scratch)
                except OSError:
                    base.errors.append(
                            ('Could not create scratch directory %s; '
                             'check that it\'s not a file and that '
                             'write permissions are active.') % scratch
                        )


    @staticmethod
    def add_args(required_parser, general_parser, output_parser, 
                    prep=False, align=False, parallel=False):
        """ Adds parameter descriptions relevant to local mode to an object
            of class argparse.ArgumentParser.

            prep: preprocess-only
            align: align-only
            parallel: add parallel-mode arguments

            No return value.
        """
        if align:
            required_parser.add_argument(
                '-i', '--input', type=str, required=True, metavar='<dir>',
                help='input directory with preprocessed reads; must be local'
            )
        if prep:
            output_parser.add_argument(
                '-o', '--output', type=str, required=False, metavar='<dir>',
                default='./rail-rna_prep',
                help='output directory; must be local or on S3'
            )
        else:
            output_parser.add_argument(
                '-o', '--output', type=str, required=False, metavar='<dir>',
                default='./rail-rna_out',
                help='output directory; must be local or on S3'
            )
        general_parser.add_argument(
            '--log', type=str, required=False, metavar='<dir>',
            default='./rail-rna_logs',
            help='directory for storing intermediate files and logs'
        )
        if not parallel:
            general_parser.add_argument(
               '-p', '--num-processes', type=int, required=False,
                metavar='<int>', default=None,
                help=('number of processes to run simultaneously (def: # cpus '
                      '- 1 if # cpus > 1; else 1)')
            )
        else:
            general_parser.add_argument(
                '--ipcontroller-json', type=str, required=False,
                metavar='<file>',
                default=None,
                help=('path to ipcontroller-client.json file '
                      '(def: IPython default for selected profile)')
            )
            general_parser.add_argument(
                '--ipython-profile', type=str, required=False, metavar='<str>',
                default=None,
                help=('connects to this IPython profile (def: default IPython '
                      'profile')
            )
            general_parser.add_argument(
                '--scratch', type=str, required=False, metavar='<dir>',
                default=None,
                help=('where to write node stdout before copying to '
                      'intermediate directory (def: securely created '
                      'temporary directory)')
            )
        general_parser.add_argument(
            '--keep-intermediates', action='store_const', const=True,
            default=False,
            help='keep intermediate files in log directory after job flow ' \
                 'is complete'
        )
        general_parser.add_argument(
            '-g', '--gzip-intermediates', action='store_const', const=True,
            default=False,
            help='compress intermediate files; slower, but saves space'
        )
        general_parser.add_argument(
           '--gzip-level', type=int, required=False, metavar='<int>',
            default=3,
            help='level of gzip compression to use for intermediates, ' \
                 'if applicable'
        )

class RailRnaElastic(object):
    """ Checks elastic-mode input parameters and relevant programs.

        Subsumes only those parameters relevant to elastic mode. Adds errors
        to base instance of RailRnaErrors.
    """
    def __init__(self, base, check_manifest=False,
        log_uri=None, ami_version='3.2.1',
        visible_to_all_users=False, tags='',
        name='Rail-RNA Job Flow',
        action_on_failure='TERMINATE_JOB_FLOW',
        hadoop_jar=None,
        master_instance_count=1, master_instance_type='c1.xlarge',
        master_instance_bid_price=None, core_instance_count=1,
        core_instance_type=None, core_instance_bid_price=None,
        task_instance_count=0, task_instance_type=None,
        task_instance_bid_price=None, ec2_key_name=None, keep_alive=False,
        termination_protected=False, no_consistent_view=False,
        intermediate_lifetime=4):

        # CLI is REQUIRED in elastic mode
        base.check_s3(reason='Rail-RNA is running in "elastic" mode')

        # Initialize possible options
        base.instance_core_counts = {
            "m1.small"    : 1,
            "m1.large"    : 2,
            "m1.xlarge"   : 4,
            "c1.medium"   : 2,
            "c1.xlarge"   : 8,
            "m2.xlarge"   : 2,
            "m2.2xlarge"  : 4,
            "m2.4xlarge"  : 8,
            "cc1.4xlarge" : 8,
            "m3.xlarge"   : 4,
            "m3.2xlarge"  : 8,
            "c3.2xlarge"  : 8,
            "c3.4xlarge"  : 16,
            "c3.8xlarge"  : 32
        }

        base.instance_mems = {
            "m1.small"    : (2*1024), #  1.7 GB
            "m1.large"    : (8*1024), #  7.5 GB
            "m1.xlarge"   : (16*1024), # 15.0 GB
            "c1.medium"   : (2*1024), #  1.7 GB
            "c1.xlarge"   : (8*1024), #  7.0 GB
            "m2.xlarge"   : (16*1024), # 17.1 GB
            "m2.2xlarge"  : (16*1024), # 34.2 GB
            "m2.4xlarge"  : (16*1024), # 68.4 GB
            "m3.xlarge"   : (15*1024),
            "m3.2xlarge"  : (30*1024),
            "cc1.4xlarge" : (16*1024), # 23.0 GB
            "c3.2xlarge" : (15*1024), # 15.0 GB
            "c3.4xlarge" : (30*1024), # 30 GB
            "c3.8xlarge" : (60*1024) # 60 GB
        }

        # From http://docs.aws.amazon.com/ElasticMapReduce/latest/
        # DeveloperGuide/TaskConfiguration_H2.html
        base.nodemanager_mems = {
            "m1.small"    : 1024,
            "m1.large"    : 3072,
            "m1.xlarge"   : 12288,
            "c1.medium"   : 1024,
            "c1.xlarge"   : 5120,
            "m2.xlarge"   : 14336,
            "m2.2xlarge"  : 30720,
            "m2.4xlarge"  : 61440,
            "cc1.4xlarge" : 20480,
            "m3.xlarge"   : 11520,
            "m3.2xlarge"  : 23040,
            "c3.2xlarge" : 11520,
            "c3.4xlarge" : 23040,
            "c3.8xlarge" : 53248
        }

        '''Not currently in use, but may become important if there are
        32- vs. 64-bit issues: base.instance_bits = {
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
            base.errors.append('Log URI (--log-uri) must be on S3, but '
                               '"{0}" was entered.'.format(log_uri))
        base.log_uri = log_uri
        base.visible_to_all_users = visible_to_all_users
        base.tags = [tag.strip() for tag in tags.split(',')]
        if len(base.tags) == 1 and base.tags[0] == '':
            base.tags = []
        base.name = name
        base.ami_version = ami_version
        # Initialize ansible for easy checks
        ansible = ab.Ansible(aws_exe=base.aws_exe, profile=base.profile)
        output_dir_url = ab.Url(base.output_dir)
        if not output_dir_url.is_s3:
            base.errors.append(('Output directory (--output) must be on S3 '
                                'when running Rail-RNA in "elastic" '
                                'mode, but {0} was entered.').format(
                                        base.output_dir
                                    ))
        if base.intermediate_dir is None:
            base.intermediate_dir = base.output_dir + '.intermediate'
        intermediate_dir_url = ab.Url(base.intermediate_dir)
        if intermediate_dir_url.is_local:
            base.errors.append(('Intermediate directory (--intermediate) '
                                'must be on HDFS or S3 when running Rail-RNA '
                                'in "elastic" mode, '
                                'but {0} was entered.').format(
                                        base.intermediate_dir
                                    ))
        elif intermediate_dir_url.is_s3:
            if not (isinstance(intermediate_lifetime, int) and
                        intermediate_lifetime != 0):
                base.errors.append(('Intermediate lifetime '
                                    '(--intermediate-lifetime) must be '
                                    '-1 or > 0, but {0} was entered.').format(
                                            intermediate_lifetime
                                        ))
            else:
                # Set up rule on S3 for deleting intermediate dir
                final_intermediate_dir = intermediate_dir_url.to_url() + '/'
                while final_intermediate_dir[-2] == '/':
                    final_intermediate_dir = final_intermediate_dir[:-1]
                ansible.s3_ansible.expire_prefix(final_intermediate_dir,
                                                    days=intermediate_lifetime)
        if ansible.s3_ansible.is_dir(base.output_dir):
            if not base.force:
                base.errors.append(('Output directory {0} exists on S3, and '
                                    '--force was not invoked to permit '
                                    'overwriting it.').format(base.output_dir))
            else:
                ansible.s3_ansible.remove_dir(base.output_dir)
        # Check manifest; download it if necessary
        manifest_url = ab.Url(base.manifest)
        if manifest_url.is_curlable \
            and 'Curl' not in base.checked_programs:
            base.curl_exe = base.check_program('curl', 'Curl', '--curl',
                                    entered_exe=base.curl_exe,
                                    reason='the manifest file is on the web')
            ansible.curl_exe = base.curl_exe
        if not ansible.exists(manifest_url.to_url()):
            base.errors.append(('Manifest file (--manifest) {0} '
                                'does not exist. Check the URL and '
                                'try again.').format(base.manifest))
        else:
            if not manifest_url.is_local:
                temp_manifest_dir = tempfile.mkdtemp()
                import atexit
                from tempdel import remove_temporary_directories
                atexit.register(remove_temporary_directories,
                                    [temp_manifest_dir])
                manifest = os.path.join(temp_manifest_dir, 'MANIFEST')
                ansible.get(base.manifest, destination=manifest)
            else:
                manifest = manifest_url.to_url()
            files_to_check = []
            base.sample_count = 0
            with open(manifest) as manifest_stream:
                for line in manifest_stream:
                    if line[0] == '#' or not line.strip(): continue
                    base.sample_count += 1
                    tokens = line.strip().split('\t')
                    check_sample_label = True
                    if len(tokens) == 5:
                        files_to_check.extend([tokens[0], tokens[2]])
                    elif len(tokens) == 3:
                        files_to_check.append(tokens[0])
                    else:
                        check_sample_label = False
                        base.errors.append(('The following line from the '
                                            'manifest file {0} '
                                            'has an invalid number of '
                                            'tokens:\n{1}'
                                            ).format(
                                                    manifest_url.to_url(),
                                                    line
                                                ))
                    if check_sample_label and tokens[-1].count('-') != 2:
                        base.errors.append(('The following line from the '
                                            'manifest file {0} '
                                            'has an invalid sample label: '
                                            '\n{1}\nA valid sample label '
                                            'takes the following form:\n'
                                            '<Group ID>-<BioRep ID>-'
                                            '<TechRep ID>'
                                            ).format(
                                                    manifest_url.to_url(),
                                                    line
                                                ))
            if files_to_check:
                if check_manifest:
                    file_count = len(files_to_check)
                    # Check files in manifest only if in preprocess job flow
                    for k, filename in enumerate(files_to_check):
                        if sys.stdout.isatty():
                            sys.stdout.write(
                                    '\r\x1b[KChecking that file %d/%d '
                                    'from manifest file exists...' % (
                                                                k+1,
                                                                file_count
                                                            )
                                )
                            sys.stdout.flush()
                        filename_url = ab.Url(filename)
                        if filename_url.is_curlable \
                            and 'Curl' not in base.checked_programs:
                            base.curl_exe = base.check_program('curl', 'Curl',
                                                '--curl',
                                                entered_exe=base.curl_exe,
                                                reason=('at least one sample '
                                                  'FASTA/FASTQ from the '
                                                  'manifest file is on '
                                                  'the web'))
                            ansible.curl_exe = base.curl_exe
                        if not ansible.exists(filename_url.to_url()):
                            base.errors.append(('The file {0} from the '
                                                'manifest file {1} does not '
                                                'exist; check the URL and try '
                                                'again.').format(
                                                        filename,
                                                        manifest_url.to_url()
                                                    ))
                    if sys.stdout.isatty():
                        sys.stdout.write(
                                '\r\x1b[KChecked all files listed in manifest '
                                'file.\n'
                            )
                        sys.stdout.flush()
            else:
                base.errors.append(('Manifest file (--manifest) {0} '
                                    'has no valid lines.').format(
                                                        manifest_url.to_url()
                                                    ))
            if not manifest_url.is_s3 and output_dir_url.is_s3:
                # Copy manifest file to S3 before job flow starts
                base.manifest = path_join(True, base.output_dir + '.manifest',
                                                'MANIFEST')
                ansible.put(manifest, base.manifest)
            if not manifest_url.is_local:
                # Clean up
                shutil.rmtree(temp_manifest_dir)

        actions_on_failure \
            = set(['TERMINATE_JOB_FLOW', 'CANCEL_AND_WAIT', 'CONTINUE',
                    'TERMINATE_CLUSTER'])

        if action_on_failure not in actions_on_failure:
            base.errors.append('Action on failure (--action-on-failure) '
                               'must be one of {{"TERMINATE_JOB_FLOW", '
                               '"CANCEL_AND_WAIT", "CONTINUE", '
                               '"TERMINATE_CLUSTER"}}, but '
                               '{0} was entered.'.format(
                                                action_on_failure
                                            ))
        base.action_on_failure = action_on_failure
        if hadoop_jar is None:
            base.hadoop_jar = _hadoop_streaming_jar
        else:
            base.hadoop_jar = hadoop_jar
        instance_type_message = ('Instance type (--instance-type) must be '
                                 'in the set {{"m1.small", "m1.large", '
                                 '"m1.xlarge", "c1.medium", "c1.xlarge", '
                                 '"m2.xlarge", "m2.2xlarge", "m2.4xlarge", '
                                 '"cc1.4xlarge"}}, but {0} was entered.')
        if master_instance_type not in base.instance_core_counts:
            base.errors.append(('Master instance type '
                               '(--master-instance-type) not valid. %s')
                                % instance_type_message.format(
                                                        master_instance_type
                                                    ))
        base.master_instance_type = master_instance_type
        if core_instance_type is None:
            base.core_instance_type = base.master_instance_type
        else:
            if core_instance_type not in base.instance_core_counts:
                base.errors.append(('Core instance type '
                                    '(--core-instance-type) not valid. %s')
                                    % instance_type_message.format(
                                                        core_instance_type
                                                    ))
            base.core_instance_type = core_instance_type
        if task_instance_type is None:
            base.task_instance_type = base.master_instance_type
        else:
            if task_instance_type not in base.instance_core_counts:
                base.errors.append(('Task instance type '
                                    '(--task-instance-type) not valid. %s')
                                    % instance_type_message.format(
                                                        task_instance_type
                                                    ))
            base.task_instance_type = task_instance_type
        if master_instance_bid_price is None:
            base.spot_master = False
        else:
            if not (master_instance_bid_price > 0):
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
            if not (core_instance_bid_price > 0):
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
            if not (task_instance_bid_price > 0):
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
                               '(--master-instance-count) must be an '
                               'integer >= 1, but {0} was entered.'.format(
                                                    master_instance_count
                                                ))
        base.master_instance_count = master_instance_count
        if not (isinstance(core_instance_count, int)
                 and core_instance_count >= 1):
            base.errors.append('Core instance count '
                               '(--core-instance-count) must be an '
                               'integer >= 1, but {0} was entered.'.format(
                                                    core_instance_count
                                                ))
        base.core_instance_count = core_instance_count
        if not (isinstance(task_instance_count, int)
                and task_instance_count >= 0):
            base.errors.append('Task instance count '
                               '(--task-instance-count) must be an '
                               'integer >= 1, but {0} was entered.'.format(
                                                    task_instance_count
                                                ))
        base.task_instance_count = task_instance_count
        if base.core_instance_count > 0:
            base.mem \
                = base.instance_mems[base.core_instance_type]
            base.nodemanager_mem \
                = base.nodemanager_mems[base.core_instance_type]
            base.max_tasks \
                = base.instance_core_counts[base.core_instance_type]
            base.total_cores = (base.core_instance_count
                * base.instance_core_counts[base.core_instance_type]
                + base.task_instance_count
                * base.instance_core_counts[base.task_instance_type])
        else:
            base.mem \
                = base.instance_mems[base.master_instance_type]
            base.nodemanager_mem \
                = base.nodemanager_mems[base.master_instance_type]
            base.max_tasks \
                = base.instance_core_counts[base.master_instance_type]
            base.total_cores \
                = base.instance_core_counts[base.master_instance_type]
        base.ec2_key_name = ec2_key_name
        base.keep_alive = keep_alive
        base.termination_protected = termination_protected
        base.no_consistent_view = no_consistent_view

    @staticmethod
    def add_args(general_parser, required_parser, output_parser, 
                    elastic_parser, align=False):
        if align:
            required_parser.add_argument(
                '-i', '--input', type=str, required=True, metavar='<s3_dir>',
                help='input directory with preprocessed read; must begin ' \
                     'with s3://'
            )
        required_parser.add_argument(
            '-o', '--output', type=str, required=True, metavar='<s3_dir>',
            help='output directory; must begin with s3://'
        )
        general_parser.add_argument(
            '--intermediate', type=str, required=False,
            metavar='<s3_dir/hdfs_dir>',
            default=None,
            help='directory for storing intermediate files; can begin with ' \
                 'hdfs:// or s3://; use S3 and set --intermediate-lifetime ' \
                 'to -1 to keep intermediates (def: output directory + ' \
                 '".intermediate")'
        )
        elastic_parser.add_argument(
            '--intermediate-lifetime', type=int, required=False,
            metavar='<int>',
            default=4,
            help='create rule for deleting intermediate files on S3 in ' \
                 'specified number of days; use -1 to keep intermediates'
        )
        elastic_parser.add_argument('--name', type=str, required=False,
            metavar='<str>',
            default='Rail-RNA Job Flow',
            help='job flow name'
        )
        elastic_parser.add_argument('--log-uri', type=str, required=False,
            metavar='<s3_dir>',
            default=None,
            help=('Hadoop log directory on S3 (def: output directory + '
                  '".logs")')
        )
        elastic_parser.add_argument('--ami-version', type=str, required=False,
            metavar='<str>',
            default='3.2.1',
            help='Amazon Machine Image to use'
        )
        elastic_parser.add_argument('--visible-to-all-users',
            action='store_const',
            const=True,
            default=False,
            help='make EC2 cluster visible to all IAM users within EMR CLI'
        )
        elastic_parser.add_argument('--action-on-failure', type=str,
            required=False,
            metavar='<choice>',
            default='TERMINATE_JOB_FLOW',
            help=('action to take if job flow fails on a given step. '
                  '<choice> is in {"TERMINATE_JOB_FLOW", "CANCEL_AND_WAIT", '
                  '"CONTINUE", "TERMINATE_CLUSTER"}')
        )
        elastic_parser.add_argument('--no-consistent-view',
            action='store_const',
            const=True,
            default=False,
            help=('do not use "consistent view," which incurs DynamoDB '
                 'charges; some intermediate data may then (very rarely) '
                 'be lost'))
        elastic_parser.add_argument('--hadoop-jar', type=str, required=False,
            metavar='<jar>',
            default=None,
            help=('Hadoop Streaming Java ARchive to use (def: AMI default)')
        )
        elastic_parser.add_argument('--master-instance-count', type=int,
            metavar='<int>',
            required=False,
            default=1,
            help=('number of master instances')
        )
        required_parser.add_argument('-c', '--core-instance-count', type=int,
            metavar='<int>',
            required=True,
            help=('number of core instances')
        )
        elastic_parser.add_argument('--task-instance-count', type=int,
            metavar='<int>',
            required=False,
            default=0,
            help=('number of task instances')
        )
        elastic_parser.add_argument('--master-instance-bid-price', type=float,
            metavar='<dec>',
            required=False,
            default=None,
            help=('bid price (dollars/hr); invoke only if master instances '
                  'should be spot')
        )
        elastic_parser.add_argument('--core-instance-bid-price', type=float,
            metavar='<dec>',
            required=False,
            default=None,
            help=('bid price (dollars/hr); invoke only if core instances '
                  'should be spot')
        )
        elastic_parser.add_argument('--task-instance-bid-price', type=float,
            metavar='<dec>',
            required=False,
            default=None,
            help=('bid price (dollars/hr); invoke only if task instances '
                  'should be spot')
        )
        elastic_parser.add_argument('--master-instance-type', type=str,
            metavar='<choice>',
            required=False,
            default='c3.2xlarge',
            help=('master instance type')
        )
        elastic_parser.add_argument('--core-instance-type', type=str,
            metavar='<choice>',
            required=False,
            default=None,
            help=('core instance type')
        )
        elastic_parser.add_argument('--task-instance-type', type=str,
            metavar='<choice>',
            required=False,
            default=None,
            help=('task instance type')
        )
        elastic_parser.add_argument('--ec2-key-name', type=str,
            metavar='<str>',
            required=False,
            default=None,
            help=('key pair name for SSHing to EC2 instances (def: '
                  'unspecified, so SSHing is not permitted)')
        )
        elastic_parser.add_argument('--keep-alive', action='store_const',
            const=True,
            default=False,
            help=('keep cluster alive after job flow completes')
        )
        elastic_parser.add_argument('--termination-protected',
            action='store_const',
            const=True,
            default=False,
            help=('protect cluster from termination in case of step failure')
        )
        elastic_parser.add_argument('--region', type=str,
            metavar='<choice>',
            required=False,
            default='us-east-1',
            help=('Amazon data center in which to run job flow. Google '
                  '"Elastic MapReduce regions" for recent list of centers ')
        )

    @staticmethod
    def hadoop_debugging_steps(base):
        return [
            {
                'ActionOnFailure' : base.action_on_failure,
                'HadoopJarStep' : {
                    'Args' : [
                        ('s3://us-east-1.elasticmapreduce/libs/'
                         'state-pusher/0.1/fetch')
                    ],
                    'Jar' : ('s3://us-east-1.elasticmapreduce/libs/'
                             'script-runner/script-runner.jar')
                },
                'Name' : 'Set up Hadoop Debugging'
            }
        ]

    @staticmethod
    def misc_steps(base):
        return [
            {
                'ActionOnFailure' : base.action_on_failure,
                'HadoopJarStep' : {
                    'Args' : [
                        '--src,s3://rail-emr/index/hg19_UCSC.tar.gz',
                        '--dest,hdfs:///index/'
                    ],
                    'Jar' : '/home/hadoop/lib/emr-s3distcp-1.0.jar'
                },
                'Name' : 'Copy Bowtie indexes from S3 to HDFS'
            }
        ]

    @staticmethod
    def bootstrap(base):
        return [
            {
                'Name' : 'Allocate swap space',
                'ScriptBootstrapAction' : {
                    'Args' : [
                        '%d' % base.mem
                    ],
                    'Path' : 's3://elasticmapreduce/bootstrap-actions/add-swap'
                }
            },
            {
                'Name' : 'Configure Hadoop',
                'ScriptBootstrapAction' : {
                    'Args' : [
                        '-c',
                        'fs.s3n.multipart.uploads.enabled=true',
                        '-y',
                        'yarn.nodemanager.pmem-check-enabled=false',
                        '-y',
                        'yarn.nodemanager.vmem-check-enabled=false',
                        '-y',
                        'yarn.nodemanager.resource.memory-mb=%d'
                        % base.nodemanager_mem,
                        '-y',
                        'yarn.scheduler.minimum-allocation-mb=%d'
                        % (base.nodemanager_mem / base.max_tasks * 8 / 10),
                        '-y',
                        'yarn.nodemanager.vmem-pmem-ratio=2.1',
                        '-y',
                        'yarn.nodemanager.container-manager.thread-count=1',
                        '-y',
                        'yarn.nodemanager.localizer.fetch.thread-count=1',
                        '-m',
                        'mapreduce.map.speculative=true',
                        '-m',
                        'mapreduce.reduce.speculative=true',
                        '-m',
                        'mapreduce.map.memory.mb=%d'
                        % (base.nodemanager_mem / base.max_tasks * 8 / 10),
                        '-m',
                        'mapreduce.reduce.memory.mb=%d'
                        % (base.nodemanager_mem / base.max_tasks * 8 / 10),
                        '-m',
                        'mapreduce.map.java.opts=-Xmx%dm'
                        % (base.nodemanager_mem / base.max_tasks * 8 / 10),
                        '-m',
                        'mapreduce.reduce.java.opts=-Xmx%dm'
                        % (base.nodemanager_mem / base.max_tasks * 8 / 10),
                        '-m',
                        'mapreduce.map.cpu.vcores=1',
                        '-m',
                        'mapreduce.reduce.cpu.vcores=1',
                        '-m',
                        'mapred.output.compress=true',
                        '-m',
                        ('mapreduce.output.fileoutputformat.compress.codec='
                         'com.hadoop.compression.lzo.LzopCodec'),
                        '-m',
                        'mapreduce.job.maps=%d' % base.total_cores,
                    ] + (['-e', 'fs.s3.consistent=true']
                            if not base.no_consistent_view else []),
                    'Path' : ('s3://elasticmapreduce/bootstrap-actions/'
                              'configure-hadoop')
                }
            }
        ]

    @staticmethod
    def instances(base):
        assert base.master_instance_count >= 1
        to_return = {
            'HadoopVersion' : '2.4.0',
            'InstanceGroups' : [
                {
                    'InstanceCount' : base.master_instance_count,
                    'InstanceRole' : 'MASTER',
                    'InstanceType': base.master_instance_type,
                    'Name' : 'Master Instance Group'
                }
            ],
            'KeepJobFlowAliveWhenNoSteps': ('true' if base.keep_alive
                                               else 'false'),
            'TerminationProtected': ('true' if base.termination_protected
                                        else 'false')
        }
        if base.master_instance_bid_price is not None:
            to_return['InstanceGroups'][0]['BidPrice'] \
                = '%0.03f' % base.master_instance_bid_price
            to_return['InstanceGroups'][0]['Market'] \
                = 'SPOT'
        else:
            to_return['InstanceGroups'][0]['Market'] \
                = 'ON_DEMAND'
        if base.core_instance_count:
            to_return['InstanceGroups'].append(
                    {
                        'InstanceCount' : base.core_instance_count,
                        'InstanceRole' : 'CORE',
                        'InstanceType': base.core_instance_type,
                        'Name' : 'Core Instance Group'
                    }
                )
            if base.core_instance_bid_price is not None:
                to_return['InstanceGroups'][1]['BidPrice'] \
                    = '%0.03f' % base.core_instance_bid_price
                to_return['InstanceGroups'][1]['Market'] \
                    = 'SPOT'
            else:
                to_return['InstanceGroups'][1]['Market'] \
                    = 'ON_DEMAND'
        if base.task_instance_count:
            to_return['InstanceGroups'].append(
                    {
                        'InstanceCount' : base.task_instance_count,
                        'InstanceRole' : 'TASK',
                        'InstanceType': base.task_instance_type,
                        'Name' : 'Task Instance Group'
                    }
                )
            if base.task_instance_bid_price is not None:
                to_return['InstanceGroups'][1]['BidPrice'] \
                    = '%0.03f' % base.task_instance_bid_price
                to_return['InstanceGroups'][1]['Market'] \
                    = 'SPOT'
            else:
                to_return['InstanceGroups'][1]['Market'] \
                    = 'ON_DEMAND'
        if base.ec2_key_name is not None:
            to_return['Ec2KeyName'] = base.ec2_key_name
        return to_return

class RailRnaPreprocess(object):
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
    def add_args(general_parser, output_parser, elastic=False):
        """ Adds parameter descriptions relevant to preprocess job flow to an
            object of class argparse.ArgumentParser.

            No return value.
        """
        if not elastic:
            output_parser.add_argument(
                '--nucleotides-per-input', type=int, required=False,
                metavar='<int>',
                default=100000000,
                help='max nucleotides from input reads to assign to each task'
            )
            output_parser.add_argument(
                '--do-not-gzip-input', action='store_const', const=True,
                default=False,
                help=('leave preprocessed input reads uncompressed')
            )
        general_parser.add_argument(
            '--do-not-check-manifest', action='store_const', const=True,
            default=False,
            help='do not check that files listed in manifest file exist'
        )

    @staticmethod
    def protosteps(base, prep_dir, push_dir, elastic=False):
        return [
            {
                'name' : 'Preprocess input reads',
                'run' : ('preprocess.py --nucs-per-file={0} {1} '
                         '--push={2} {3}').format(
                                                    base.nucleotides_per_input,
                                                    '--gzip-output' if
                                                    base.gzip_input else '',
                                                    ab.Url(push_dir).to_url(
                                                            caps=True
                                                        ),
                                                    '--stdout' if elastic
                                                    else ''
                                                ),
                'inputs' : [base.manifest],
                'no_input_prefix' : True,
                'output' : push_dir if elastic else prep_dir,
                'no_output_prefix' : True,
                'inputformat' : (
                       'org.apache.hadoop.mapred.lib.NLineInputFormat'
                    ),
                'taskx' : 0,
                'index_output' : True,
                'direct_copy' : True
            },
        ]

    @staticmethod
    def bootstrap():
        return [
            {
                'Name' : 'Install PyPy',
                'ScriptBootstrapAction' : {
                    'Args' : [
                        ('s3://rail-emr/bin/'
                         'pypy-2.2.1-linux_x86_64-portable.tar.bz2')
                    ],  
                    'Path' : 's3://rail-emr/bootstrap/install-pypy.sh'
                }
            },
            {
                'Name' : 'Install Rail-RNA',
                'ScriptBootstrapAction' : {
                    'Args' : [
                        's3://rail-emr/bin/rail-rna-%s.tar.gz' 
                        % version_number,
                        '/mnt'
                    ],
                    'Path' : 's3://rail-emr/bootstrap/install-rail.sh'
                }
            }
        ]

class RailRnaAlign(object):
    """ Sets parameters relevant to just the "align" job flow. """
    def __init__(self, base, input_dir=None, elastic=False,
        bowtie1_exe=None, bowtie1_idx='genome', bowtie1_build_exe=None,
        bowtie2_exe=None, bowtie2_build_exe=None, bowtie2_idx='genome',
        bowtie2_args='', samtools_exe=None, bedgraphtobigwig_exe=None,
        genome_partition_length=5000, max_readlet_size=25,
        readlet_config_size=32, min_readlet_size=15, readlet_interval=4,
        cap_size_multiplier=1.2, max_intron_size=500000, min_intron_size=10,
        min_exon_size=9, search_filter='none',
        motif_search_window_size=1000, max_gaps_mismatches=3,
        motif_radius=5, genome_bowtie1_args='-v 0 -a -m 80',
        transcriptome_bowtie2_args='-k 30', count_multiplier=6, tie_margin=6,
        normalize_percentile=0.75, very_replicable=False,
        do_not_output_bam_by_chr=False, output_sam=False,
        bam_basename='alignments', bed_basename='', assembly='hg19',
        s3_ansible=None):
        if not elastic:
            '''Programs and Bowtie indices should be checked only in local
            mode.'''
            base.bowtie1_exe = base.check_program('bowtie', 'Bowtie 1',
                                '--bowtie1', entered_exe=bowtie1_exe,
                                is_exe=is_exe, which=which)
            bowtie1_version_command = [base.bowtie1_exe, '--version']
            try:
                base.bowtie1_version = subprocess.check_output(
                        bowtie1_version_command
                    ).split('\n', 1)[0].split(' ')[-1]
            except Exception as e:
                base.errors.append(('Error "{0}" encountered attempting to '
                                    'execute "{1}".').format(
                                                        e.message,
                                                        ' '.join(
                                                        bowtie1_version_command
                                                       )
                                                    ))
            base.bowtie1_build_exe = base.check_program('bowtie-build',
                                            'Bowtie 1 Build',
                                            '--bowtie1-build',
                                            entered_exe=bowtie1_build_exe,
                                            is_exe=is_exe,
                                            which=which)
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
            base.bowtie1_idx = bowtie1_idx
            base.bowtie2_exe = base.check_program('bowtie2', 'Bowtie 2',
                                '--bowtie2', entered_exe=bowtie2_exe,
                                is_exe=is_exe, which=which)
            bowtie2_version_command = [base.bowtie2_exe, '--version']
            try:
                base.bowtie2_version = subprocess.check_output(
                        bowtie2_version_command
                    ).split('\n', 1)[0].split(' ')[-1]
            except Exception as e:
                base.errors.append(('Error "{0}" encountered attempting to '
                                    'execute "{1}".').format(
                                                        e.message,
                                                        ' '.join(
                                                        bowtie2_version_command
                                                       )
                                                    ))
            base.bowtie2_build_exe = base.check_program('bowtie2-build',
                                            'Bowtie 2 Build',
                                            '--bowtie2-build',
                                            entered_exe=bowtie2_build_exe,
                                            is_exe=is_exe,
                                            which=which)
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
            base.bowtie2_idx = bowtie2_idx
            base.samtools_exe = base.check_program('samtools', 'SAMTools',
                                '--samtools', entered_exe=samtools_exe,
                                is_exe=is_exe, which=which)
            try:
                samtools_process = subprocess.Popen(
                        [base.samtools_exe], stderr=subprocess.PIPE
                    )
            except Exception as e:
                base.errors.append(('Error "{0}" encountered attempting to '
                                    'execute "{1}".').format(
                                                        e.message,
                                                        base.samtools_exe
                                                    ))
            if base.errors:
                # Output any errors before detect message is determined
                raise RuntimeError(
                        '\n'.join(
                                ['%d) %s' % (i+1, error) for i, error
                                    in enumerate(base.errors)]
                            ) if len(base.errors) > 1 else base.errors[0]
                    )
            base.samtools_version = '<unknown>'
            for line in samtools_process.stderr:
                if 'Version:' in line:
                    base.samtools_version = line.rpartition(' ')[-1].strip()
                    if base.samtools_version[-1] == ')':
                        base.samtools_version = base.samtools_version[:-1]
            base.detect_message =('Detected Bowtie 1 v{0}, Bowtie 2 v{1}, '
                                  'and SAMTools v{2}.').format(
                                               base.bowtie1_version,
                                               base.bowtie2_version,
                                               base.samtools_version
                                            )
            base.bedgraphtobigwig_exe = base.check_program('bedGraphToBigWig', 
                                    'BedGraphToBigWig', '--bedgraphtobigwig',
                                    entered_exe=bedgraphtobigwig_exe,
                                    is_exe=is_exe, which=which)
            # Check input dir
            if input_dir is not None:
                if not os.path.exists(input_dir):
                    base_errors.append(('Input directory (--input) '
                                        '"{0}" does not exist').format(
                                                            input_dir
                                                        ))
                else:
                    base.input_dir = input_dir
        else:
            # Elastic mode; check S3 for genome if necessary
            assert s3_ansible is not None
            if assembly == 'hg19':
                base.index_archive = 's3://rail-emr/index/hg19_UCSC.tar.gz'
            else:
                if not Url(assembly).is_s3:
                    base.errors.append(('Bowtie index archive must be on S3'
                                        ' in "elastic" mode, but '
                                        '"{0}" was entered.').format(assembly))
                elif not s3_ansible.exists(assembly):
                    base.errors.append('Bowtie index archive was not found '
                                       'on S3 at "{0}".'.format(assembly))
                else:
                    base.index_archive = assembly
            if input_dir is not None:
                if not ab.Url(input_dir).is_s3:
                    base.errors.append(('Input directory must be on S3, but '
                                        '"{0}" was entered.').format(
                                                                input_dir
                                                            ))
                elif not s3_ansible.is_dir(input_dir):
                    base.errors.append(('Input directory "{0}" was not found '
                                        'on S3.').format(input_dir))
                else:
                    base.input_dir = input_dir
            # Set up elastic params
            base.bowtie1_idx = _elastic_bowtie1_idx
            base.bowtie2_idx = _elastic_bowtie2_idx
            base.bedgraphtobigwig_exe = _elastic_bedgraphtobigwig_exe
            base.samtools_exe = _elastic_samtools_exe
            base.bowtie1_exe = _elastic_bowtie1_exe
            base.bowtie2_exe = _elastic_bowtie2_exe
            base.bowtie1_build_exe = _elastic_bowtie1_build_exe
            base.bowtie2_build_exe = _elastic_bowtie2_build_exe

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
        if not (isinstance(readlet_config_size, int) and readlet_config_size
                >= max_readlet_size):
            base.errors.append('Readlet config size (--readlet-config-size) '
                               'must be an integer >= maximum readlet size '
                               '(--max-readlet-size) = '
                               '{0}, but {1} was entered.'.format(
                                                    base.max_readlet_size,
                                                    readlet_config_size
                                                ))
        base.readlet_config_size = readlet_config_size
        if not (isinstance(readlet_interval, int) and readlet_interval
                > 0):
            base.errors.append('Readlet interval (--readlet-interval) '
                               'must be an integer > 0, '
                               'but {0} was entered.'.format(
                                                    readlet_interval
                                                ))
        base.readlet_interval = readlet_interval
        if not (cap_size_multiplier > 1):
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
        if search_filter == 'none':
            base.search_filter = 0
        elif search_filter == 'mild':
            base.search_filter = base.min_exon_size
        elif search_filter == 'strict':
            try:
                base.search_filter = int(base.min_exon_size * 1.5)
            except:
                pass
        elif not (
                (isinstance(search_filter, int)
                    and search_filter >= 0) or
                search_filter in ['none', 'mild', 'strict']
            ):
            base.errors.append('Search filter (--search-filter) '
                               'must be an integer >= 0 or one of {"none", '
                               '"mild", "strict"}, but {0} was '
                               'entered.'.format(search_filter))
        else:
            base.search_filter = search_filter
        if not (isinstance(motif_search_window_size, int) and 
                    motif_search_window_size >= 0):
            base.errors.append('Motif search window size '
                               '(--motif-search-window-size) must be an '
                               'integer >= 0, but {0} was entered.'.format(
                                                    motif_search_window_size
                                                ))
        base.motif_search_window_size = motif_search_window_size
        if max_gaps_mismatches is not None and not (
                isinstance(max_gaps_mismatches, int) and 
                max_gaps_mismatches >= 0
            ):
            base.errors.append('Max gaps and mismatches '
                               '(--max-gaps-mismatches) must be an '
                               'integer >= 0, but {0} was entered.'.format(
                                                    max_gaps_mismatches
                                                ))
        base.max_gaps_mismatches = max_gaps_mismatches
        if not (isinstance(motif_radius, int) and
                    motif_radius >= 0):
            base.errors.append('Motif radius (--motif-radius) must be an '
                               'integer >= 0, but {0} was entered.'.format(
                                                    motif_radius
                                                ))
        base.motif_radius = motif_radius
        base.genome_bowtie1_args = genome_bowtie1_args
        base.transcriptome_bowtie2_args = transcriptome_bowtie2_args
        if not (0 <= normalize_percentile <= 1):
            base.errors.append('Normalization percentile '
                               '(--normalize-percentile) must on the '
                               'interval [0, 1], but {0} was entered'.format(
                                                    normalize_percentile
                                                ))
        base.normalize_percentile = normalize_percentile
        if not (isinstance(tie_margin, int) and
                    tie_margin >= 0):
            base.errors.append('Tie margin (--tie-margin) must be an '
                               'integer >= 0, but {0} was entered.'.format(
                                                    tie_margin
                                                ))
        base.tie_margin = tie_margin
        if not (isinstance(count_multiplier, int) and
                    count_multiplier >= 0):
            base.errors.append('Count multiplier (--count-multiplier) must '
                               'be an integer >= 0, but '
                               '{0} was entered.'.format(
                                                    count_multiplier
                                                ))
        base.count_multiplier = count_multiplier
        base.do_not_output_bam_by_chr = do_not_output_bam_by_chr
        base.very_replicable = very_replicable
        base.output_sam = output_sam
        base.bam_basename = bam_basename
        base.bed_basename = bed_basename

    @staticmethod
    def add_args(required_parser, exec_parser, output_parser, algo_parser, 
                    elastic=False):
        """ usage: argparse.SUPPRESS if advanced options should be suppressed;
                else None
        """
        if not elastic:
            exec_parser.add_argument(
                '--bowtie1', type=str, required=False,
                metavar='<exe>',
                default=None,
                help=('path to Bowtie 1 executable (def: bowtie)')
            )
            exec_parser.add_argument(
                '--bowtie1-build', type=str, required=False,
                metavar='<exe>',
                default=None,
                help=('path to Bowtie 1 Build executable (def: bowtie-build)')
            )
            required_parser.add_argument(
                '-1', '--bowtie1-idx', type=str, required=True,
                metavar='<idx>',
                help='path to Bowtie 1 index; include basename'
            )
            exec_parser.add_argument(
                '--bowtie2', type=str, required=False,
                metavar='<exe>',
                default=None,
                help=('path to Bowtie 2 executable (def: bowtie2)')
            )
            exec_parser.add_argument(
                '--bowtie2-build', type=str, required=False,
                metavar='<exe>',
                default=None,
                help=('path to Bowtie 2 Build executable (def: bowtie2-build)')
            )
            required_parser.add_argument(
                '-2', '--bowtie2-idx', type=str, required=True,
                metavar='<idx>',
                help='path to Bowtie 2 index; include basename'
            )
            exec_parser.add_argument(
                '--samtools', type=str, required=False,
                metavar='<exe>',
                default=None,
                help=('path to SAMTools executable (def: samtools)')
            )
            exec_parser.add_argument(
                '--bedgraphtobigwig', type=str, required=False,
                metavar='<exe>',
                default=None,
                help=('path to BedGraphToBigWig executable '
                      '(def: bedGraphToBigWig)')
            )
        else:
            required_parser.add_argument(
                '-a', '--assembly', type=str, required=True,
                metavar='<choice/tgz>',
                help=('assembly to use for alignment. <choice> can be in '
                      '{"hg19"}. otherwise, specify path to tar.gz Rail '
                      'archive on S3')
            )
        algo_parser.add_argument(
                '--bowtie2-args', type=str, required=False,
                default='',
                metavar='<str>',
                help=('arguments to pass to Bowtie 2, which is always run in '
                      '"--local" mode (def: Bowtie 2 defaults)')
            )
        algo_parser.add_argument(
            '--genome-partition-length', type=int, required=False,
            metavar='<int>',
            default=5000,
            help=('smallest unit of genome addressable by single task when '
                  'computing coverage')
        )
        algo_parser.add_argument(
            '--max-readlet-size', type=int, required=False,
            metavar='<int>',
            default=25,
            help='max size of read segment to align when searching for introns'
        )
        algo_parser.add_argument(
            '--readlet-config-size', type=int, required=False,
            metavar='<int>',
            default=35,
            help=('max number of exonic bases spanned by a path enumerated in '
                  'intron DAG')
        )
        algo_parser.add_argument(
            '--min-readlet-size', type=int, required=False,
            metavar='<int>',
            default=15,
            help='min size of read segment to align when searching for introns'
        )
        algo_parser.add_argument(
            '--readlet-interval', type=int, required=False,
            metavar='<int>',
            default=4,
            help=('distance between start positions of successive overlapping '
                  'read segments to align when searching for introns')
        )
        algo_parser.add_argument(
            '--cap-size-multiplier', type=float, required=False,
            default=1.1,
            help=SUPPRESS
        )
        algo_parser.add_argument(
            '--max-intron-size', type=int, required=False,
            metavar='<int>',
            default=500000,
            help=('filter out introns spanning more than <int> bases')
        )
        algo_parser.add_argument(
            '--min-intron-size', type=int, required=False,
            metavar='<int>',
            default=10,
            help=('filter out introns spanning fewer than <int> bases')
        )
        algo_parser.add_argument(
            '--min-exon-size', type=int, required=False,
            metavar='<int>',
            default=9,
            help=('try to be sensitive to exons that span at least <int> '
                  'bases')
        )
        algo_parser.add_argument(
            '--search-filter', type=str, required=False,
            metavar='<choice/int>',
            default='none',
            help=('filter out reads searched for introns that fall below '
                  'threshold <int> for initially detected anchor length; '
                  'or select <choice> from {"strict", "mild", "none"}')
        )
        algo_parser.add_argument(
            '--motif-search-window-size', type=int, required=False,
            default=1000,
            help=SUPPRESS
        )
        algo_parser.add_argument(
            '--max-gaps-mismatches', type=int, required=False,
            default=None,
            help=SUPPRESS
        )
        algo_parser.add_argument(
            '--motif-radius', type=int, required=False,
            default=5,
            help=SUPPRESS
        )
        algo_parser.add_argument(
            '--genome-bowtie1-args', type=str, required=False,
            default='-v 0 -a -m 30',
            help=SUPPRESS
        )
        algo_parser.add_argument(
            '--transcriptome-bowtie2-args', type=str, required=False,
            default='-k 30',
            help=SUPPRESS
        )
        algo_parser.add_argument(
            '--count-multiplier', type=int, required=False,
            default=15,
            help=SUPPRESS
        )
        algo_parser.add_argument(
            '--normalize-percentile', type=float, required=False,
            metavar='<dec>',
            default=0.75,
            help=('percentile to use when computing normalization factors for '
                  'sample coverages')
        )
        algo_parser.add_argument(
            '--tie-margin', type=int, required=False,
            metavar='<int>',
            default=0,
            help=('allowed score difference per 100 bases among ties in '
                  'max score. For example, 150 and 144 are tied alignment '
                  'scores for a 100-bp read when --tie-margin is 6')
        )
        algo_parser.add_argument(
            '--very-replicable', action='store_const', const=True,
            default=False,
            help=('ensures that results are exactly replicable across '
                  'modes, core counts, and cluster sizes; otherwise, '
                  'replicability is nearly exact unless same configuration '
                  'is used')
        )
        output_parser.add_argument(
            '--do-not-output-bam-by-chr', action='store_const', const=True,
            default=False,
            help=('place all of a sample\'s alignments in one file rather '
                  'than dividing them up by chromosome')
        )
        output_parser.add_argument(
            '--output-sam', action='store_const', const=True,
            default=False,
            help='output SAM instead of BAM'
        )
        output_parser.add_argument(
            '--bam-basename', type=str, required=False,
            metavar='<str>',
            default='alignments',
            help='basename for BAM output'
        )
        output_parser.add_argument(
            '--bed-basename', type=str, required=False,
            metavar='<str>',
            default='',
            help='basename for BED output (def: *empty*)'
        )

    @staticmethod
    def protosteps(base, input_dir, elastic=False):
        manifest = ('/mnt/MANIFEST' if elastic else base.manifest)
        verbose = ('--verbose' if base.verbose else '')
        keep_alive = ('--keep-alive' if elastic else '')
        return [  
            {
                'name' : 'Align reads to genome',
                'run' : ('align_reads.py --bowtie-idx={0} --bowtie2-idx={1} '
                         '--bowtie2-exe={2} '
                         '--exon-differentials --partition-length={3} '
                         '--min-exon-size={4} '
                         '--manifest={5} {6} {7} {8} -- {9}').format(
                                                        base.bowtie1_idx,
                                                        base.bowtie2_idx,
                                                        base.bowtie2_exe,
                                                base.genome_partition_length,
                                                base.search_filter,
                                                        manifest,
                                                        verbose,
                                                        keep_alive,
                                                        '--ignore-first-token'
                                                        if elastic else '',
                                                        base.bowtie2_args),
                'inputs' : [input_dir],
                'no_input_prefix' : True,
                'output' : 'align_reads',
                'taskx' : 0,
                'multiple_outputs' : True,
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size / 4)
                    ]
            },
            {
                'name' : 'Segment unique read sequences into readlets',
                'run' : ('readletize.py --max-readlet-size={0} '
                         '--readlet-interval={1} '
                         '--capping-multiplier={2}').format(
                                    base.max_readlet_size,
                                    base.readlet_interval,
                                    base.cap_size_multiplier
                                ),
                'inputs' : [path_join(elastic, 'align_reads', 'readletize')],
                'output' : 'readletize',
                'taskx' : 3,
                'part' : 'k1,1',
                'keys' : 1,
                'multiple_outputs' : True,
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % _base_combine_split_size
                    ]
            },
            {
                'name' : 'Align unique readlets to genome',
                'run' : ('align_readlets.py --bowtie-idx={0} '
                         '--bowtie-exe={1} {2} {3} '
                         '-- -t --sam-nohead --startverbose {4}').format(
                                                    base.bowtie1_idx,
                                                    base.bowtie1_exe,
                                                    verbose,
                                                    keep_alive,
                                                    base.genome_bowtie1_args,
                                                ),
                'inputs' : [path_join(elastic, 'readletize', 'readletized')],
                'output' : 'align_readlets',
                'taskx' : 3,
                'part' : 'k1,1',
                'keys' : 1,
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size*4),
                        'mapreduce.reduce.memory.mb=%d'
                        % ((base.nodemanager_mem / base.max_tasks * 8 / 10 * 2)
                            if hasattr(base, 'nodemanager_mem') else 0),
                    ],
                'direct_copy' : True
            },
            {
                'name' : 'Search for introns using readlet alignments',
                'run' : ('intron_search.py --bowtie-idx={0} '
                         '--partition-length={1} --max-intron-size={2} '
                         '--min-intron-size={3} --min-exon-size={4} '
                         '--search-window-size={5} {6} '
                         '--motif-radius={7} {8}').format(
                                                base.bowtie1_idx,
                                                base.genome_partition_length,
                                                base.max_intron_size,
                                                base.min_intron_size,
                                                base.min_exon_size,
                                                base.motif_search_window_size,
                                                ('--max-gaps-mismatches %d' %
                                                 base.max_gaps_mismatches)
                                                if base.max_gaps_mismatches
                                                is not None else '',
                                                base.motif_radius,
                                                verbose
                                            ),
                'inputs' : ['align_readlets'],
                'output' : 'intron_search',
                'taskx' : 3,
                'part' : 'k1,1',
                'keys' : 1,
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size*4)
                    ],
                'direct_copy' : True
            },
            {
                'name' : 'Enumerate possible intron cooccurrences on readlets',
                'run' : ('intron_config.py '
                         '--readlet-size={0} {1}').format(
                                                    base.readlet_config_size,
                                                    verbose
                                                ),
                'inputs' : ['intron_search'],
                'output' : 'intron_config',
                'taskx' : 1,
                'part' : 'k1,2',
                'keys' : 4,
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size)
                    ],
                'direct_copy' : True
            },
            {
                'name' : 'Get transcriptome elements for realignment',
                'run' : ('intron_fasta.py --bowtie-idx={0} {1}').format(
                                                        base.bowtie1_idx,
                                                        verbose
                                                    ),
                'inputs' : ['intron_config'],
                'output' : 'intron_fasta',
                'taskx' : 1,
                'part' : 'k1,4',
                'keys' : 4,
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size)
                    ],
                'direct_copy' : True
            },
            {
                'name' : 'Build index of transcriptome elements',
                'run' : ('intron_index.py --bowtie-build-exe={0} '
                         '--out={1} {2}').format(base.bowtie2_build_exe,
                                                 ab.Url(
                                                    path_join(elastic,
                                                        base.output_dir,
                                                        'transcript_index')
                                                    ).to_url(caps=True),
                                                 keep_alive),
                'inputs' : ['intron_fasta'],
                'output' : 'intron_index',
                'taskx' : None,
                'part' : 'k1,1',
                'keys' : 1,
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size)
                    ],
                'direct_copy' : True
            },
            {
                'name' : 'Finalize intron cooccurrences on reads',
                'run' : ('cointron_enum.py --bowtie2-idx={0} '
                         '--bowtie2-exe={1} {2} {3} -- {4}').format(
                                            'intron/intron'
                                            if elastic else
                                            path_join(elastic,
                                                base.output_dir,
                                                'transcript_index',
                                                'intron'),
                                            base.bowtie2_exe,
                                            verbose,
                                            keep_alive,
                                            base.transcriptome_bowtie2_args
                                        ),
                'inputs' : [path_join(elastic, 'readletize', 'unique')],
                'output' : 'cointron_enum',
                'taskx' : max(3, base.sample_count / 30),
                'archives' : ab.Url(path_join(elastic,
                                    base.output_dir,
                                    'transcript_index',
                                    'intron.tar.gz#intron')).to_native_url(),
                'part' : 'k1,1',
                'keys' : 1,
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size*2)
                    ],
                'direct_copy' : True
            },
            {
                'name' : 'Get transcriptome elements for read realignment',
                'run' : ('cointron_fasta.py --bowtie-idx={0} {1}').format(
                                                        base.bowtie1_idx,
                                                        verbose
                                                    ),
                'inputs' : ['cointron_enum'],
                'output' : 'cointron_fasta',
                'taskx' : max(8, base.sample_count / 30),
                'part' : 'k1,4',
                'keys' : 7,
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size)
                    ],
                'direct_copy' : True
            },
            {
                'name' : 'Align reads to transcriptome elements',
                'run' : ('realign_reads.py --bowtie2-exe={0} '
                         '--count-multiplier {1} {2} {3} {4} -- {5}').format(
                                        base.bowtie2_exe,
                                        base.count_multiplier,
                                        verbose,
                                        keep_alive,
                                        '--replicable'
                                        if base.very_replicable
                                        else '',
                                        base.bowtie2_args
                                    ),
                'inputs' : [path_join(elastic, 'align_reads', 'unmapped'),
                            'cointron_fasta'],
                'output' : 'realign_reads',
                # Ensure that a single reducer isn't assigned too much fasta
                'taskx' : max(6, base.sample_count / 20),
                'part' : 'k1,1',
                'keys' : 1,
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size*3)
                    ],
                'direct_copy' : True
            },
            {
                'name' : 'Collect and compare read alignments',
                'run' : ('compare_alignments.py --bowtie-idx={0} '
                         '--partition-length={1} --exon-differentials '
                         '--tie-margin {2} --manifest={3} {4} -- {5}').format(
                                        base.bowtie1_idx,
                                        base.genome_partition_length,
                                        base.tie_margin,
                                        manifest,
                                        verbose,
                                        base.bowtie2_args
                                    ),
                'inputs' : [path_join(elastic, 'align_reads', 'postponed_sam'),
                            'realign_reads'],
                'output' : 'compare_alignments',
                # Ensure that a single reducer isn't assigned too much fasta
                'taskx' : max(4, base.sample_count / 20),
                'part' : 'k1,1',
                'keys' : 1,
                'multiple_outputs' : True,
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size*3)
                    ]
            },
            {
                'name' : 'Associate spliced alignments with intron coverages',
                'run' : 'intron_coverage.py --bowtie-idx {0}'.format(
                                                        base.bowtie1_idx
                                                    ),
                'inputs' : [path_join(elastic, 'compare_alignments',
                                                    'intron_bed'),
                            path_join(elastic, 'compare_alignments',
                                               'sam_intron_ties')],
                'output' : 'intron_coverage',
                'taskx' : 1,
                'part' : 'k1,6',
                'keys' : 7,
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size)
                    ],
                'direct_copy' : True
            },
            {
                'name' : 'Finalize primary alignments of spliced reads',
                'run' : ('break_ties.py --exon-differentials '
                            '--bowtie-idx {0} --partition-length {1} '
                            '--manifest {2} -- {3}').format(
                                    base.bowtie1_idx,
                                    base.genome_partition_length,
                                    manifest,
                                    base.bowtie2_args
                                ),
                'inputs' : ['intron_coverage',
                            path_join(elastic, 'compare_alignments',
                                               'sam_clip_ties')],
                'output' : 'break_ties',
                'taskx' : 1,
                'part' : 'k1,1',
                'keys' : 1,
                'multiple_outputs' : True,
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size)
                    ]
            },
            {
                'name' : 'Merge exon differentials at same genomic positions',
                'run' : 'sum.py {0}'.format(
                                        keep_alive
                                    ),
                'inputs' : [path_join(elastic, 'align_reads', 'exon_diff'),
                            path_join(elastic, 'compare_alignments',
                                               'exon_diff'),
                            path_join(elastic, 'break_ties', 'exon_diff')],
                'output' : 'collapse',
                'taskx' : max(4, base.sample_count / 20),
                'part' : 'k1,3',
                'keys' : 3,
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size*2)
                    ],
                'direct_copy' : True
            },
            {
                'name' : 'Compile sample coverages from exon differentials',
                'run' : ('coverage_pre.py --bowtie-idx={0} '
                         '--partition-stats').format(base.bowtie1_idx),
                'inputs' : ['collapse'],
                'output' : 'precoverage',
                'taskx' : max(4, base.sample_count / 20),
                'part' : 'k1,2',
                'keys' : 3,
                'multiple_outputs' : True,
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size*2)
                    ]
            },
            {
                'name' : 'Write bigwigs with exome coverage by sample',
                'run' : ('coverage.py --bowtie-idx={0} --percentile={1} '
                         '--out={2} --bigwig-exe={3} '
                         '--manifest={4} {5}').format(base.bowtie1_idx,
                                                     base.normalize_percentile,
                                                     ab.Url(
                                                        path_join(elastic,
                                                        base.output_dir,
                                                        'coverage_bigwigs')
                                                     ).to_url(caps=True),
                                                     base.bedgraphtobigwig_exe,
                                                     manifest,
                                                     verbose),
                'inputs' : [path_join(elastic, 'precoverage', 'coverage')],
                'output' : 'coverage',
                'taskx' : 1,
                'part' : 'k1,1',
                'keys' : 3,
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size)
                    ],
                'direct_copy' : True
            },
            {
                'name' : 'Write normalization factors for sample coverages',
                'run' : 'coverage_post.py --out={0} --manifest={1}'.format(
                                                        ab.Url(
                                                            path_join(elastic,
                                                            base.output_dir,
                                                    'normalization_factors')
                                                        ).to_url(caps=True),
                                                        manifest
                                                    ),
                'inputs' : ['coverage'],
                'output' : 'coverage_post',
                'taskx' : None,
                'part' : 'k1,1',
                'keys' : 2,
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size)
                    ],
                'direct_copy' : True
            },
            {
                'name' : 'Aggregate intron and indel results by sample',
                'run' : 'bed_pre.py',
                'inputs' : [path_join(elastic, 'compare_alignments',
                                               'indel_bed'),
                            path_join(elastic, 'align_reads', 'indel_bed'),
                            path_join(elastic, 'break_ties', 'indel_bed'),
                            path_join(elastic, 'compare_alignments',
                                               'intron_bed'),
                            path_join(elastic, 'break_ties', 'intron_bed')],
                'output' : 'prebed',
                'taskx' : max(4, base.sample_count / 20),
                'part' : 'k1,6',
                'keys' : 6,
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size*2)
                    ],
                'direct_copy' : True
            },
            {
                'name' : 'Write beds with intron and indel results by sample',
                'run' : ('bed.py --bowtie-idx={0} --out={1} '
                         '--manifest={2} --bed-basename={3}').format(
                                                        base.bowtie1_idx,
                                                        ab.Url(
                                                            path_join(elastic,
                                                            base.output_dir,
                                                        'introns_and_indels')
                                                         ).to_url(caps=True),
                                                        manifest,
                                                        base.bed_basename
                                                    ),
                'inputs' : ['prebed'],
                'output' : 'bed',
                'taskx' : 1,
                'part' : 'k1,2',
                'keys' : 4,
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size)
                    ],
                'direct_copy' : True
            },
            {
                'name' : 'Write bams with alignments by sample',
                'run' : ('bam.py --out={0} --bowtie-idx={1} '
                         '--samtools-exe={2} --bam-basename={3} '
                         '--manifest={4} {5} {6}').format(
                                        ab.Url(
                                            path_join(elastic,
                                            base.output_dir, 'alignments')
                                        ).to_url(caps=True),
                                        base.bowtie1_idx,
                                        base.samtools_exe,
                                        base.bam_basename,
                                        manifest,
                                        keep_alive,
                                        '--output-by-chromosome'
                                        if not base.do_not_output_bam_by_chr
                                        else ''
                                    ),
                'inputs' : [path_join(elastic, 'align_reads', 'sam'),
                            path_join(elastic, 'compare_alignments', 'sam'),
                            path_join(elastic, 'break_ties', 'sam')],
                'output' : 'bam',
                'taskx' : 1,
                'part' : ('k1,1' if base.do_not_output_bam_by_chr else 'k1,2'),
                'keys' : 3,
                'extra_args' : [
                        'mapreduce.reduce.shuffle.input.buffer.percent=0.4',
                        'mapreduce.reduce.shuffle.merge.percent=0.4',
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size)
                    ],
                'direct_copy' : True
            }]

    @staticmethod
    def bootstrap(base):
        return [
            {
                'Name' : 'Install PyPy',
                'ScriptBootstrapAction' : {
                    'Args' : [],
                    'Path' : 's3://rail-emr/bootstrap/install-pypy.sh'
                }
            },
            {
                'Name' : 'Install Bowtie 1',
                'ScriptBootstrapAction' : {
                    'Args' : [],
                    'Path' : 's3://rail-emr/bootstrap/install-bowtie.sh'
                }
            },
            {
                'Name' : 'Install Bowtie 2',
                'ScriptBootstrapAction' : {
                    'Args' : [],
                    'Path' : 's3://rail-emr/bootstrap/install-bowtie2.sh'
                }
            },
            {
                'Name' : 'Install bedGraphToBigWig',
                'ScriptBootstrapAction' : {
                    'Args' : [
                        '/mnt/bin'
                    ],
                    'Path' : 's3://rail-emr/bootstrap/install-kenttools.sh'
                }
            },
            {
                'Name' : 'Install SAMTools',
                'ScriptBootstrapAction' : {
                    'Args' : [],
                    'Path' : 's3://rail-emr/bootstrap/install-samtools.sh'
                }
            },
            {
                'Name' : 'Install Rail-RNA',
                'ScriptBootstrapAction' : {
                    'Args' : [
                        's3://rail-emr/bin/rail-rna-%s.tar.gz'
                        % version_number,
                        '/mnt'
                    ],
                    'Path' : 's3://rail-emr/bootstrap/install-rail.sh'
                }
            },
            {
                'Name' : 'Transfer Bowtie indexes to nodes',
                'ScriptBootstrapAction' : {
                    'Args' : [
                        '/mnt',
                        base.index_archive
                    ],
                    'Path' : 's3://rail-emr/bootstrap/install-index.sh'
                }
            },
            {
                'Name' : 'Transfer manifest file to nodes',
                'ScriptBootstrapAction' : {
                    'Args' : [
                        base.manifest,
                        '/mnt',
                        'MANIFEST'
                    ],
                    'Path' : 's3://rail-emr/bootstrap/s3cmd_s3.sh'
                }
            }
        ]

class RailRnaLocalPreprocessJson(object):
    """ Constructs JSON for local mode + preprocess job flow. """
    def __init__(self, manifest, output_dir, intermediate_dir='./intermediate',
        force=False, aws_exe=None, profile='default', region='us-east-1',
        verbose=False, nucleotides_per_input=8000000, gzip_input=True,
        num_processes=1, gzip_intermediates=False, gzip_level=3,
        keep_intermediates=False, check_manifest=True):
        base = RailRnaErrors(manifest, output_dir, 
            intermediate_dir=intermediate_dir,
            force=force, aws_exe=aws_exe, profile=profile,
            region=region, verbose=verbose)
        RailRnaLocal(base, check_manifest=check_manifest,
            num_processes=num_processes, gzip_intermediates=gzip_intermediates,
            gzip_level=gzip_level, keep_intermediates=keep_intermediates)
        RailRnaPreprocess(base, nucleotides_per_input=nucleotides_per_input,
            gzip_input=gzip_input)
        if base.errors:
            raise RuntimeError(
                    '\n'.join(
                            ['%d) %s' % (i+1, error) for i, error
                                in enumerate(base.errors)]
                        ) if len(base.errors) > 1 else base.errors[0]
                )
        self._json_serial = {}
        step_dir = os.path.join(base_path, 'rna', 'steps')
        self._json_serial['Steps'] = steps(RailRnaPreprocess.protosteps(base,
                os.path.join(base.intermediate_dir, 'preprocess'),
                base.output_dir, elastic=False),
                '', '', step_dir,
                base.num_processes,
                base.intermediate_dir, unix=False
            )
        self.base = base
    
    @property
    def json_serial(self):
        return self._json_serial

class RailRnaParallelPreprocessJson(object):
    """ Constructs JSON for parallel mode + preprocess job flow. """
    def __init__(self, manifest, output_dir, intermediate_dir='./intermediate',
        force=False, aws_exe=None, profile='default', region='us-east-1',
        verbose=False, nucleotides_per_input=8000000, gzip_input=True,
        num_processes=1, gzip_intermediates=False, gzip_level=3,
        ipython_profile=None, ipcontroller_json=None, scratch=None,
        keep_intermediates=False, check_manifest=True):
        rc = ipython_client(ipython_profile=ipython_profile,
                                ipcontroller_json=ipcontroller_json)
        ready_engines(rc)
        base = RailRnaErrors(manifest, output_dir, 
            intermediate_dir=intermediate_dir,
            force=force, aws_exe=aws_exe, profile=profile,
            region=region, verbose=verbose)
        engine_bases = [RailRnaErrors(manifest, output_dir, 
            intermediate_dir=intermediate_dir,
            force=force, aws_exe=aws_exe, profile=profile,
            region=region, verbose=verbose) for i in rc.ids]
        RailRnaLocal(base, check_manifest=check_manifest,
            num_processes=num_processes, gzip_intermediates=gzip_intermediates,
            gzip_level=gzip_level, keep_intermediates=keep_intermediates,
            local=False, parallel=False)
        asyncresults = []
        for i in rc.ids:
            asyncresults.append(
                    rc[i].apply_async(
                        RailRnaLocal.__init__, engine_bases[i],
                        check_manifest=check_manifest,
                        num_processes=num_processes,
                        gzip_intermediates=gzip_intermediates,
                        gzip_level=gzip_level,
                        keep_intermediates=keep_intermediates,
                        local=False, parallel=True, ansible=ab.Ansible()
                    )
                )
        while any([not asyncresult.ready() for asyncresult in asyncresults]):
            time.sleep(1e-1)
        asyncexceptions = defaultdict(set)
        for asyncresult in asyncresults:
            try:
                asyncdict = asyncresult.get_dict()
            except Exception as e:
                asyncexceptions[format_exc()].add(
                        asyncresult.metadata['engine_id']
                    )
        if asyncexceptions:
            runtimeerror_message = []
            for exc in asyncexceptions:
                runtimeerror_message.extend(
                        ['Engine(s) %s report(s) the following exception.'
                            % list(asyncexceptions[exc]),
                         exc]
                     )
            raise RuntimeError('\n'.join(runtimeerror_message))
        asyncresults = []
        if base.check_curl_on_engines:
            for i in rc.ids:
                asyncresults.append(
                        rc[i].apply_async(
                            engine_bases[i].check_program, 'curl', 'Curl',
                            '--curl', entered_exe=base.curl_exe,
                            reason=base.check_curl_on_engines,
                            is_exe=is_exe,
                            which=which
                        )
                    )
        if base.check_s3_on_engines:
            for i in rc.ids:
                asyncresults.append(
                        rc[i].apply_async(
                            engine_bases[i].check_s3,
                            reason=base.check_curl_on_engines,
                            is_exe=is_exe,
                            which=which
                        )
                    )
        if asyncresults:
            asyncexceptions = defaultdict(set)
            for asyncresult in asyncresults:
                try:
                    asyncdict = asyncresult.get_dict()
                except Exception as e:
                    asyncexceptions[format_exc()].add(
                            asyncresult.metadata['engine_id']
                        )
            if asyncexceptions:
                runtimeerror_message = []
                for exc in asyncexceptions:
                    runtimeerror_message.extend(
                            ['Engines %s report the following exception.'
                                % asyncexceptions[exc],
                             exc]
                         )
                raise RuntimeError('\n'.join(runtimeerror_message))
        RailRnaPreprocess(base,
            nucleotides_per_input=nucleotides_per_input, gzip_input=gzip_input)
        if base.errors:
            raise RuntimeError(
                    '\n'.join(
                            ['%d) %s' % (i+1, error) for i, error
                                in enumerate(base.errors)]
                        ) if len(base.errors) > 1 else base.errors[0]
                )
        self._json_serial = {}
        step_dir = os.path.join(base_path, 'rna', 'steps')
        self._json_serial['Steps'] = steps(RailRnaPreprocess.protosteps(base,
                os.path.join(base.intermediate_dir, 'preprocess'),
                base.output_dir, elastic=False),
                '', '', step_dir,
                base.num_processes,
                base.intermediate_dir, unix=False
            )
        self.base = base
    
    @property
    def json_serial(self):
        return self._json_serial

class RailRnaElasticPreprocessJson(object):
    """ Constructs JSON for elastic mode + preprocess job flow. """
    def __init__(self, manifest, output_dir, intermediate_dir='./intermediate',
        force=False, aws_exe=None, profile='default', region='us-east-1',
        verbose=False, nucleotides_per_input=8000000, gzip_input=True,
        log_uri=None, ami_version='3.2.1',
        visible_to_all_users=False, tags='',
        name='Rail-RNA Job Flow',
        action_on_failure='TERMINATE_JOB_FLOW',
        hadoop_jar=None,
        master_instance_count=1, master_instance_type='c1.xlarge',
        master_instance_bid_price=None, core_instance_count=1,
        core_instance_type=None, core_instance_bid_price=None,
        task_instance_count=0, task_instance_type=None,
        task_instance_bid_price=None, ec2_key_name=None, keep_alive=False,
        termination_protected=False, no_consistent_view=False,
        check_manifest=True, intermediate_lifetime=4):
        base = RailRnaErrors(manifest, output_dir, 
            intermediate_dir=intermediate_dir,
            force=force, aws_exe=aws_exe, profile=profile,
            region=region, verbose=verbose)
        RailRnaElastic(base, check_manifest=check_manifest,
            log_uri=log_uri, ami_version=ami_version,
            visible_to_all_users=visible_to_all_users, tags=tags,
            name=name,
            action_on_failure=action_on_failure,
            hadoop_jar=hadoop_jar, master_instance_count=master_instance_count,
            master_instance_type=master_instance_type,
            master_instance_bid_price=master_instance_bid_price,
            core_instance_count=core_instance_count,
            core_instance_type=core_instance_type,
            core_instance_bid_price=core_instance_bid_price,
            task_instance_count=task_instance_count,
            task_instance_type=task_instance_type,
            task_instance_bid_price=task_instance_bid_price,
            ec2_key_name=ec2_key_name, keep_alive=keep_alive,
            termination_protected=termination_protected,
            no_consistent_view=no_consistent_view,
            intermediate_lifetime=intermediate_lifetime)
        RailRnaPreprocess(base, nucleotides_per_input=nucleotides_per_input,
            gzip_input=gzip_input)
        if base.errors:
            raise RuntimeError(
                    '\n'.join(
                            ['%d) %s' % (i+1, error) for i, error
                                in enumerate(base.errors)]
                        ) if len(base.errors) > 1 else base.errors[0]
                )
        self._json_serial = {}
        if base.core_instance_count > 0:
            reducer_count = base.core_instance_count \
                * base.instance_core_counts[base.core_instance_type]
        else:
            reducer_count = base.master_instance_count \
                * base.instance_core_counts[base.core_instance_type]
        self._json_serial['Steps'] \
            = RailRnaElastic.hadoop_debugging_steps(base) + steps(
                    RailRnaPreprocess.protosteps(base,
                        path_join(True, base.intermediate_dir, 'preprocess'),
                        base.output_dir, elastic=True),
                    base.action_on_failure,
                    base.hadoop_jar, '/mnt/src/rna/steps',
                    reducer_count, base.intermediate_dir, unix=True,
                    no_consistent_view=base.no_consistent_view
                )
        self._json_serial['AmiVersion'] = base.ami_version
        if base.log_uri is not None:
            self._json_serial['LogUri'] = base.log_uri
        else:
            self._json_serial['LogUri'] = base.output_dir + '.logs'
        self._json_serial['Name'] = base.name
        self._json_serial['NewSupportedProducts'] = []
        self._json_serial['Tags'] = base.tags
        self._json_serial['VisibleToAllUsers'] = (
                'true' if base.visible_to_all_users else 'false'
            )
        self._json_serial['Instances'] = RailRnaElastic.instances(base)
        self._json_serial['BootstrapActions'] \
            = RailRnaPreprocess.bootstrap() \
            + RailRnaElastic.bootstrap(base)
        self.base = base
    
    @property
    def json_serial(self):
        return self._json_serial

class RailRnaLocalAlignJson(object):
    """ Constructs JSON for local mode + align job flow. """
    def __init__(self, manifest, output_dir, input_dir,
        intermediate_dir='./intermediate',
        force=False, aws_exe=None, profile='default', region='us-east-1',
        verbose=False, bowtie1_exe=None,
        bowtie1_idx='genome', bowtie1_build_exe=None, bowtie2_exe=None,
        bowtie2_build_exe=None, bowtie2_idx='genome',
        bowtie2_args='', samtools_exe=None, bedgraphtobigwig_exe=None,
        genome_partition_length=5000, max_readlet_size=25,
        readlet_config_size=32, min_readlet_size=15, readlet_interval=4,
        cap_size_multiplier=1.2, max_intron_size=500000, min_intron_size=10,
        min_exon_size=9, search_filter='none',
        motif_search_window_size=1000, max_gaps_mismatches=3,
        motif_radius=5, genome_bowtie1_args='-v 0 -a -m 80',
        transcriptome_bowtie2_args='-k 30', count_multiplier=6, 
        tie_margin=6, very_replicable=False, normalize_percentile=0.75,
        do_not_output_bam_by_chr=False,
        output_sam=False, bam_basename='alignments',
        bed_basename='', num_processes=1, gzip_intermediates=False,
        gzip_level=3, keep_intermediates=False):
        base = RailRnaErrors(manifest, output_dir, 
            intermediate_dir=intermediate_dir,
            force=force, aws_exe=aws_exe, profile=profile,
            region=region, verbose=verbose)
        RailRnaLocal(base, check_manifest=False, num_processes=num_processes,
            gzip_intermediates=gzip_intermediates, gzip_level=gzip_level,
            keep_intermediates=keep_intermediates)
        RailRnaAlign(base, input_dir=input_dir,
            elastic=False, bowtie1_exe=bowtie1_exe,
            bowtie1_idx=bowtie1_idx, bowtie1_build_exe=bowtie1_build_exe,
            bowtie2_exe=bowtie2_exe, bowtie2_build_exe=bowtie2_build_exe,
            bowtie2_idx=bowtie2_idx, bowtie2_args=bowtie2_args,
            samtools_exe=samtools_exe,
            bedgraphtobigwig_exe=bedgraphtobigwig_exe,
            genome_partition_length=genome_partition_length,
            max_readlet_size=max_readlet_size,
            readlet_config_size=readlet_config_size,
            min_readlet_size=min_readlet_size,
            readlet_interval=readlet_interval,
            cap_size_multiplier=cap_size_multiplier,
            max_intron_size=max_intron_size,
            min_intron_size=min_intron_size, min_exon_size=min_exon_size,
            search_filter=search_filter,
            motif_search_window_size=motif_search_window_size,
            max_gaps_mismatches=max_gaps_mismatches,
            motif_radius=motif_radius,
            genome_bowtie1_args=genome_bowtie1_args,
            count_multiplier=count_multiplier,
            transcriptome_bowtie2_args=transcriptome_bowtie2_args,
            tie_margin=tie_margin,
            normalize_percentile=normalize_percentile,
            very_replicable=very_replicable,
            do_not_output_bam_by_chr=do_not_output_bam_by_chr,
            output_sam=output_sam, bam_basename=bam_basename,
            bed_basename=bed_basename)
        if base.errors:
            raise RuntimeError(
                    '\n'.join(
                            ['%d) %s' % (i+1, error) for i, error
                                in enumerate(base.errors)]
                        ) if len(base.errors) > 1 else base.errors[0]
                )
        print >>sys.stderr, base.detect_message
        if not sys.stderr.isatty():
            # So the user sees it too
            print base.detect_message
        self._json_serial = {}
        step_dir = os.path.join(base_path, 'rna', 'steps')
        self._json_serial['Steps'] = steps(RailRnaAlign.protosteps(base,
                base.input_dir, elastic=False), '', '', step_dir,
                base.num_processes, base.intermediate_dir, unix=False
            )
        self.base = base

    @property
    def json_serial(self):
        return self._json_serial

class RailRnaParallelAlignJson(object):
    """ Constructs JSON for local mode + align job flow. """
    def __init__(self, manifest, output_dir, input_dir,
        intermediate_dir='./intermediate',
        force=False, aws_exe=None, profile='default', region='us-east-1',
        verbose=False, bowtie1_exe=None,
        bowtie1_idx='genome', bowtie1_build_exe=None, bowtie2_exe=None,
        bowtie2_build_exe=None, bowtie2_idx='genome',
        bowtie2_args='', samtools_exe=None, bedgraphtobigwig_exe=None,
        genome_partition_length=5000, max_readlet_size=25,
        readlet_config_size=32, min_readlet_size=15, readlet_interval=4,
        cap_size_multiplier=1.2, max_intron_size=500000, min_intron_size=10,
        min_exon_size=9, search_filter='none',
        motif_search_window_size=1000, max_gaps_mismatches=3,
        motif_radius=5, genome_bowtie1_args='-v 0 -a -m 80',
        transcriptome_bowtie2_args='-k 30', count_multiplier=6, 
        tie_margin=6, very_replicable=False, normalize_percentile=0.75,
        do_not_output_bam_by_chr=False,
        output_sam=False, bam_basename='alignments',
        bed_basename='', num_processes=1,
        ipython_profile=None, ipcontroller_json=None, scratch=None,
        gzip_intermediates=False,
        gzip_level=3, keep_intermediates=False):
        rc = ipython_client(ipython_profile=ipython_profile,
                                ipcontroller_json=ipcontroller_json)
        ready_engines(rc)
        engine_bases = [RailRnaErrors(manifest, output_dir, 
            intermediate_dir=intermediate_dir,
            force=force, aws_exe=aws_exe, profile=profile,
            region=region, verbose=verbose) for i in rc.ids]
        RailRnaLocal(base, check_manifest=False, num_processes=num_processes,
            gzip_intermediates=gzip_intermediates, gzip_level=gzip_level,
            keep_intermediates=keep_intermediates, local=False, parallel=False)
        asyncresults = []
        for i in rc.ids:
            asyncresults.append(
                    rc[i].apply_async(
                        RailRnaLocal.__init__, engine_bases[i],
                        check_manifest=check_manifest,
                        num_processes=num_processes,
                        gzip_intermediates=gzip_intermediates,
                        gzip_level=gzip_level,
                        keep_intermediates=keep_intermediates,
                        local=False, parallel=True, ansible=ab.Ansible()
                    )
                )
        while any([not asyncresult.ready() for asyncresult in asyncresults]):
            time.sleep(1e-1)
        asyncexceptions = defaultdict(set)
        for asyncresult in asyncresults:
            try:
                asyncdict = asyncresult.get_dict()
            except Exception as e:
                asyncexceptions[format_exc()].add(
                        asyncresult.metadata['engine_id']
                    )
        if asyncexceptions:
            runtimeerror_message = []
            for exc in asyncexceptions:
                runtimeerror_message.extend(
                        ['Engine(s) %s report(s) the following exception.'
                            % list(asyncexceptions[exc]),
                         exc]
                     )
            raise RuntimeError('\n'.join(runtimeerror_message))
        asyncresults = []
        if base.check_curl_on_engines:
            for i in rc.ids:
                asyncresults.append(
                        rc[i].apply_async(
                            engine_bases[i].check_program, 'curl', 'Curl',
                            '--curl', entered_exe=base.curl_exe,
                            reason=base.check_curl_on_engines,
                            is_exe=is_exe,
                            which=which
                        )
                    )
        if base.check_s3_on_engines:
            for i in rc.ids:
                asyncresults.append(
                        rc[i].apply_async(
                            engine_bases[i].check_s3,
                            reason=base.check_curl_on_engines,
                            is_exe=is_exe,
                            which=which
                        )
                    )
        if asyncresults:
            asyncexceptions = defaultdict(set)
            for asyncresult in asyncresults:
                try:
                    asyncdict = asyncresult.get_dict()
                except Exception as e:
                    asyncexceptions[format_exc()].add(
                            asyncresult.metadata['engine_id']
                        )
            if asyncexceptions:
                runtimeerror_message = []
                for exc in asyncexceptions:
                    runtimeerror_message.extend(
                            ['Engine(s) %s report(s) the following exception.'
                                % list(asyncexceptions[exc]),
                             exc]
                         )
                raise RuntimeError('\n'.join(runtimeerror_message))
        RailRnaAlign(base, input_dir=input_dir,
            elastic=False, bowtie1_exe=bowtie1_exe,
            bowtie1_idx=bowtie1_idx, bowtie1_build_exe=bowtie1_build_exe,
            bowtie2_exe=bowtie2_exe, bowtie2_build_exe=bowtie2_build_exe,
            bowtie2_idx=bowtie2_idx, bowtie2_args=bowtie2_args,
            samtools_exe=samtools_exe,
            bedgraphtobigwig_exe=bedgraphtobigwig_exe,
            genome_partition_length=genome_partition_length,
            max_readlet_size=max_readlet_size,
            readlet_config_size=readlet_config_size,
            min_readlet_size=min_readlet_size,
            readlet_interval=readlet_interval,
            cap_size_multiplier=cap_size_multiplier,
            max_intron_size=max_intron_size,
            min_intron_size=min_intron_size, min_exon_size=min_exon_size,
            search_filter=search_filter,
            motif_search_window_size=motif_search_window_size,
            max_gaps_mismatches=max_gaps_mismatches,
            motif_radius=motif_radius,
            genome_bowtie1_args=genome_bowtie1_args,
            count_multiplier=count_multiplier,
            transcriptome_bowtie2_args=transcriptome_bowtie2_args,
            tie_margin=tie_margin,
            normalize_percentile=normalize_percentile,
            very_replicable=very_replicable,
            do_not_output_bam_by_chr=do_not_output_bam_by_chr,
            output_sam=output_sam, bam_basename=bam_basename,
            bed_basename=bed_basename)
        if base.errors:
            raise RuntimeError(
                    '\n'.join(
                            ['%d) %s' % (i+1, error) for i, error
                                in enumerate(base.errors)]
                        ) if len(base.errors) > 1 else base.errors[0]
                )
        print >>sys.stderr, base.detect_message
        if not sys.stderr.isatty():
            # So the user sees it too
            print base.detect_message
        self._json_serial = {}
        step_dir = os.path.join(base_path, 'rna', 'steps')
        self._json_serial['Steps'] = steps(RailRnaAlign.protosteps(base,
                base.input_dir, elastic=False), '', '', step_dir,
                base.num_processes, base.intermediate_dir, unix=False
            )
        self.base = base

    @property
    def json_serial(self):
        return self._json_serial

class RailRnaElasticAlignJson(object):
    """ Constructs JSON for elastic mode + align job flow. """
    def __init__(self, manifest, output_dir, input_dir, 
        intermediate_dir='./intermediate',
        force=False, aws_exe=None, profile='default', region='us-east-1',
        verbose=False, bowtie1_exe=None, bowtie1_idx='genome',
        bowtie1_build_exe=None, bowtie2_exe=None,
        bowtie2_build_exe=None, bowtie2_idx='genome',
        bowtie2_args='', samtools_exe=None, bedgraphtobigwig_exe=None,
        genome_partition_length=5000, max_readlet_size=25,
        readlet_config_size=32, min_readlet_size=15, readlet_interval=4,
        cap_size_multiplier=1.2, max_intron_size=500000, min_intron_size=10,
        min_exon_size=9, search_filter='none',
        motif_search_window_size=1000, max_gaps_mismatches=3,
        motif_radius=5, genome_bowtie1_args='-v 0 -a -m 80',
        transcriptome_bowtie2_args='-k 30', count_multiplier=6,
        tie_margin=6, very_replicable=False,
        normalize_percentile=0.75, do_not_output_bam_by_chr=False,
        output_sam=False, bam_basename='alignments',
        bed_basename='', log_uri=None, ami_version='3.2.1',
        visible_to_all_users=False, tags='',
        name='Rail-RNA Job Flow',
        action_on_failure='TERMINATE_JOB_FLOW',
        hadoop_jar=None,
        master_instance_count=1, master_instance_type='c1.xlarge',
        master_instance_bid_price=None, core_instance_count=1,
        core_instance_type=None, core_instance_bid_price=None,
        task_instance_count=0, task_instance_type=None,
        task_instance_bid_price=None, ec2_key_name=None, keep_alive=False,
        termination_protected=False, no_consistent_view=False,
        intermediate_lifetime=4):
        base = RailRnaErrors(manifest, output_dir, 
            intermediate_dir=intermediate_dir,
            force=force, aws_exe=aws_exe, profile=profile,
            region=region, verbose=verbose)
        RailRnaElastic(base, check_manifest=False,
            log_uri=log_uri, ami_version=ami_version,
            visible_to_all_users=visible_to_all_users, tags=tags,
            name=name, action_on_failure=action_on_failure,
            hadoop_jar=hadoop_jar, master_instance_count=master_instance_count,
            master_instance_type=master_instance_type,
            master_instance_bid_price=master_instance_bid_price,
            core_instance_count=core_instance_count,
            core_instance_type=core_instance_type,
            core_instance_bid_price=core_instance_bid_price,
            task_instance_count=task_instance_count,
            task_instance_type=task_instance_type,
            task_instance_bid_price=task_instance_bid_price,
            ec2_key_name=ec2_key_name, keep_alive=keep_alive,
            termination_protected=termination_protected,
            no_consistent_view=no_consistent_view,
            intermediate_lifetime=intermediate_lifetime)
        RailRnaAlign(base, input_dir=input_dir,
            elastic=True, bowtie1_exe=bowtie1_exe,
            bowtie1_idx=bowtie1_idx, bowtie1_build_exe=bowtie1_build_exe,
            bowtie2_exe=bowtie2_exe, bowtie2_build_exe=bowtie2_build_exe,
            bowtie2_idx=bowtie2_idx, bowtie2_args=bowtie2_args,
            samtools_exe=samtools_exe,
            bedgraphtobigwig_exe=bedgraphtobigwig_exe,
            genome_partition_length=genome_partition_length,
            max_readlet_size=max_readlet_size,
            readlet_config_size=readlet_config_size,
            min_readlet_size=min_readlet_size,
            readlet_interval=readlet_interval,
            cap_size_multiplier=cap_size_multiplier,
            max_intron_size=max_intron_size,
            min_intron_size=min_intron_size, min_exon_size=min_exon_size,
            search_filter=search_filter,
            motif_search_window_size=motif_search_window_size,
            max_gaps_mismatches=max_gaps_mismatches,
            motif_radius=motif_radius,
            genome_bowtie1_args=genome_bowtie1_args,
            count_multiplier=count_multiplier,
            transcriptome_bowtie2_args=transcriptome_bowtie2_args,
            tie_margin=tie_margin,
            normalize_percentile=normalize_percentile,
            very_replicable=very_replicable,
            do_not_output_bam_by_chr=do_not_output_bam_by_chr,
            output_sam=output_sam, bam_basename=bam_basename,
            bed_basename=bed_basename,
            s3_ansible=ab.S3Ansible(aws_exe=base.aws_exe,
                                        profile=base.profile))
        if base.errors:
            raise RuntimeError(
                    '\n'.join(
                            ['%d) %s' % (i+1, error) for i, error
                                in enumerate(base.errors)]
                        ) if len(base.errors) > 1 else base.errors[0]
                )
        self._json_serial = {}
        if base.core_instance_count > 0:
            reducer_count = base.core_instance_count \
                * base.instance_core_counts[base.core_instance_type]
        else:
            reducer_count = base.master_instance_count \
                * base.instance_core_counts[base.core_instance_type]
        self._json_serial['Steps'] \
            = RailRnaElastic.hadoop_debugging_steps(base) + \
                steps(
                    RailRnaAlign.protosteps(base, base.input_dir,
                                                        elastic=True),
                    base.action_on_failure,
                    base.hadoop_jar, '/mnt/src/rna/steps',
                    reducer_count, base.intermediate_dir, unix=True,
                    no_consistent_view=base.no_consistent_view
                )
        self._json_serial['AmiVersion'] = base.ami_version
        if base.log_uri is not None:
            self._json_serial['LogUri'] = base.log_uri
        else:
            self._json_serial['LogUri'] = base.output_dir + '.logs'
        self._json_serial['Name'] = base.name
        self._json_serial['NewSupportedProducts'] = []
        self._json_serial['Tags'] = base.tags
        self._json_serial['VisibleToAllUsers'] = (
                'true' if base.visible_to_all_users else 'false'
            )
        self._json_serial['Instances'] = RailRnaElastic.instances(base)
        self._json_serial['BootstrapActions'] \
            = RailRnaAlign.bootstrap(base) \
            + RailRnaElastic.bootstrap(base)
        self.base = base
    
    @property
    def json_serial(self):
        return self._json_serial

class RailRnaLocalAllJson(object):
    """ Constructs JSON for local mode + preprocess+align job flow. """
    def __init__(self, manifest, output_dir, intermediate_dir='./intermediate',
        force=False, aws_exe=None, profile='default', region='us-east-1',
        verbose=False, nucleotides_per_input=8000000, gzip_input=True,
        bowtie1_exe=None, bowtie1_idx='genome', bowtie1_build_exe=None,
        bowtie2_exe=None, bowtie2_build_exe=None, bowtie2_idx='genome',
        bowtie2_args='', samtools_exe=None, bedgraphtobigwig_exe=None,
        genome_partition_length=5000, max_readlet_size=25,
        readlet_config_size=32, min_readlet_size=15, readlet_interval=4,
        cap_size_multiplier=1.2, max_intron_size=500000, min_intron_size=10,
        min_exon_size=9, search_filter='none',
        motif_search_window_size=1000, max_gaps_mismatches=3,
        motif_radius=5, genome_bowtie1_args='-v 0 -a -m 80',
        count_multiplier=6, transcriptome_bowtie2_args='-k 30', tie_margin=6,
        very_replicable=False, normalize_percentile=0.75,
        do_not_output_bam_by_chr=False,
        output_sam=False, bam_basename='alignments', bed_basename='',
        num_processes=1, gzip_intermediates=False, gzip_level=3,
        keep_intermediates=False, check_manifest=True):
        base = RailRnaErrors(manifest, output_dir, 
            intermediate_dir=intermediate_dir,
            force=force, aws_exe=aws_exe, profile=profile,
            region=region, verbose=verbose)
        RailRnaPreprocess(base, nucleotides_per_input=nucleotides_per_input,
            gzip_input=gzip_input)
        RailRnaLocal(base, check_manifest=check_manifest,
            num_processes=num_processes, gzip_intermediates=gzip_intermediates,
            gzip_level=gzip_level, keep_intermediates=keep_intermediates)
        RailRnaAlign(base, bowtie1_exe=bowtie1_exe,
            bowtie1_idx=bowtie1_idx, bowtie1_build_exe=bowtie1_build_exe,
            bowtie2_exe=bowtie2_exe, bowtie2_build_exe=bowtie2_build_exe,
            bowtie2_idx=bowtie2_idx, bowtie2_args=bowtie2_args,
            samtools_exe=samtools_exe,
            bedgraphtobigwig_exe=bedgraphtobigwig_exe,
            genome_partition_length=genome_partition_length,
            max_readlet_size=max_readlet_size,
            readlet_config_size=readlet_config_size,
            min_readlet_size=min_readlet_size,
            readlet_interval=readlet_interval,
            cap_size_multiplier=cap_size_multiplier,
            max_intron_size=max_intron_size,
            min_intron_size=min_intron_size, min_exon_size=min_exon_size,
            search_filter=search_filter,
            motif_search_window_size=motif_search_window_size,
            max_gaps_mismatches=max_gaps_mismatches,
            motif_radius=motif_radius,
            genome_bowtie1_args=genome_bowtie1_args,
            transcriptome_bowtie2_args=transcriptome_bowtie2_args,
            count_multiplier=count_multiplier,
            tie_margin=tie_margin,
            normalize_percentile=normalize_percentile,
            very_replicable=very_replicable,
            do_not_output_bam_by_chr=do_not_output_bam_by_chr,
            output_sam=output_sam, bam_basename=bam_basename,
            bed_basename=bed_basename)
        if base.errors:
            raise RuntimeError(
                    '\n'.join(
                            ['%d) %s' % (i+1, error) for i, error
                                in enumerate(base.errors)]
                        ) if len(base.errors) > 1 else base.errors[0]
                )
        print >>sys.stderr, base.detect_message
        if not sys.stderr.isatty():
            # So the user sees it too
            print base.detect_message
        self._json_serial = {}
        step_dir = os.path.join(base_path, 'rna', 'steps')
        prep_dir = path_join(False, base.intermediate_dir,
                                        'preprocess')
        push_dir = path_join(False, base.intermediate_dir,
                                        'preprocess', 'push')
        self._json_serial['Steps'] = \
            steps(RailRnaPreprocess.protosteps(base,
                prep_dir, push_dir, elastic=False), '', '', step_dir,
                base.num_processes, base.intermediate_dir, unix=False
            ) + \
            steps(RailRnaAlign.protosteps(base,
                push_dir, elastic=False), '', '', step_dir,
                base.num_processes, base.intermediate_dir, unix=False
            )
        self.base = base

    @property
    def json_serial(self):
        return self._json_serial

class RailRnaParallelAllJson(object):
    """ Constructs JSON for local mode + preprocess+align job flow. """
    def __init__(self, manifest, output_dir, intermediate_dir='./intermediate',
        force=False, aws_exe=None, profile='default', region='us-east-1',
        verbose=False, nucleotides_per_input=8000000, gzip_input=True,
        bowtie1_exe=None, bowtie1_idx='genome', bowtie1_build_exe=None,
        bowtie2_exe=None, bowtie2_build_exe=None, bowtie2_idx='genome',
        bowtie2_args='', samtools_exe=None, bedgraphtobigwig_exe=None,
        genome_partition_length=5000, max_readlet_size=25,
        readlet_config_size=32, min_readlet_size=15, readlet_interval=4,
        cap_size_multiplier=1.2, max_intron_size=500000, min_intron_size=10,
        min_exon_size=9, search_filter='none',
        motif_search_window_size=1000, max_gaps_mismatches=3,
        motif_radius=5, genome_bowtie1_args='-v 0 -a -m 80',
        count_multiplier=6, transcriptome_bowtie2_args='-k 30', tie_margin=6,
        very_replicable=False, normalize_percentile=0.75,
        do_not_output_bam_by_chr=False,
        output_sam=False, bam_basename='alignments', bed_basename='',
        num_processes=1, gzip_intermediates=False, gzip_level=3,
        ipython_profile=None, ipcontroller_json=None, scratch=None,
        keep_intermediates=False, check_manifest=True):
        rc = ipython_client(ipython_profile=ipython_profile,
                                ipcontroller_json=ipcontroller_json)
        ready_engines(rc)
        engine_bases = [RailRnaErrors(manifest, output_dir, 
            intermediate_dir=intermediate_dir,
            force=force, aws_exe=aws_exe, profile=profile,
            region=region, verbose=verbose) for i in rc.ids]
        base = RailRnaErrors(manifest, output_dir, 
            intermediate_dir=intermediate_dir,
            force=force, aws_exe=aws_exe, profile=profile,
            region=region, verbose=verbose)
        from IPython.parallel import require
        RailRnaLocal = require(globals()['RailRnaLocal'], rna_config)
        RailRnaLocal(base, check_manifest=check_manifest,
            num_processes=num_processes, gzip_intermediates=gzip_intermediates,
            gzip_level=gzip_level, keep_intermediates=keep_intermediates,
            parallel=False, local=False)
        asyncresults = []
        for i in rc.ids:
            asyncresults.append(
                    rc[i].apply_async(
                        RailRnaLocal, engine_bases[i],
                        check_manifest=check_manifest,
                        num_processes=num_processes,
                        gzip_intermediates=gzip_intermediates,
                        gzip_level=gzip_level,
                        keep_intermediates=keep_intermediates,
                        local=False, parallel=True, ansible=ab.Ansible()
                    )
                )
        while any([not asyncresult.ready() for asyncresult in asyncresults]):
            time.sleep(1e-1)
        asyncexceptions = defaultdict(set)
        for asyncresult in asyncresults:
            try:
                asyncdict = asyncresult.get_dict()
            except Exception as e:
                asyncexceptions[format_exc()].add(
                        asyncresult.metadata['engine_id']
                    )
        if asyncexceptions:
            runtimeerror_message = []
            for exc in asyncexceptions:
                runtimeerror_message.extend(
                        ['Engine(s) %s report(s) the following exception.'
                            % list(asyncexceptions[exc]),
                         exc]
                     )
            raise RuntimeError('\n'.join(runtimeerror_message))
        asyncresults = []
        if base.check_curl_on_engines:
            for i in rc.ids:
                asyncresults.append(
                        rc[i].apply_async(
                            engine_bases[i].check_program, 'curl', 'Curl',
                            '--curl', entered_exe=base.curl_exe,
                            reason=base.check_curl_on_engines,
                            is_exe=is_exe,
                            which=which
                        )
                    )
        if base.check_s3_on_engines:
            for i in rc.ids:
                asyncresults.append(
                        rc[i].apply_async(
                            engine_bases[i].check_s3,
                            reason=base.check_curl_on_engines,
                            is_exe=is_exe,
                            which=which
                        )
                    )
        if asyncresults:
            asyncexceptions = defaultdict(set)
            for asyncresult in asyncresults:
                try:
                    asyncdict = asyncresult.get_dict()
                except Exception as e:
                    asyncexceptions[format_exc()].add(
                            asyncresult.metadata['engine_id']
                        )
            if asyncexceptions:
                runtimeerror_message = []
                for exc in asyncexceptions:
                    runtimeerror_message.extend(
                            ['Engine(s) %s report(s) the following exception.'
                                % list(asyncexceptions[exc]),
                             exc]
                         )
                raise RuntimeError('\n'.join(runtimeerror_message))
        RailRnaPreprocess(base, nucleotides_per_input=nucleotides_per_input,
            gzip_input=gzip_input)
        RailRnaAlign(base, bowtie1_exe=bowtie1_exe,
            bowtie1_idx=bowtie1_idx, bowtie1_build_exe=bowtie1_build_exe,
            bowtie2_exe=bowtie2_exe, bowtie2_build_exe=bowtie2_build_exe,
            bowtie2_idx=bowtie2_idx, bowtie2_args=bowtie2_args,
            samtools_exe=samtools_exe,
            bedgraphtobigwig_exe=bedgraphtobigwig_exe,
            genome_partition_length=genome_partition_length,
            max_readlet_size=max_readlet_size,
            readlet_config_size=readlet_config_size,
            min_readlet_size=min_readlet_size,
            readlet_interval=readlet_interval,
            cap_size_multiplier=cap_size_multiplier,
            max_intron_size=max_intron_size,
            min_intron_size=min_intron_size,
            search_filter=search_filter,
            min_exon_size=min_exon_size,
            motif_search_window_size=motif_search_window_size,
            max_gaps_mismatches=max_gaps_mismatches,
            motif_radius=motif_radius,
            genome_bowtie1_args=genome_bowtie1_args,
            transcriptome_bowtie2_args=transcriptome_bowtie2_args,
            count_multiplier=count_multiplier,
            tie_margin=tie_margin,
            normalize_percentile=normalize_percentile,
            very_replicable=very_replicable,
            do_not_output_bam_by_chr=do_not_output_bam_by_chr,
            output_sam=output_sam, bam_basename=bam_basename,
            bed_basename=bed_basename)
        if base.errors:
            raise RuntimeError(
                    '\n'.join(
                            ['%d) %s' % (i+1, error) for i, error
                                in enumerate(base.errors)]
                        ) if len(base.errors) > 1 else base.errors[0]
                )
        print >>sys.stderr, base.detect_message
        if not sys.stderr.isatty():
            # So the user sees it too
            print base.detect_message
        self._json_serial = {}
        step_dir = os.path.join(base_path, 'rna', 'steps')
        prep_dir = path_join(False, base.intermediate_dir,
                                        'preprocess')
        push_dir = path_join(False, base.intermediate_dir,
                                        'preprocess', 'push')
        self._json_serial['Steps'] = \
            steps(RailRnaPreprocess.protosteps(base,
                prep_dir, push_dir, elastic=False), '', '', step_dir,
                base.num_processes, base.intermediate_dir, unix=False
            ) + \
            steps(RailRnaAlign.protosteps(base,
                push_dir, elastic=False), '', '', step_dir,
                base.num_processes, base.intermediate_dir, unix=False
            )
        self.base = base

    @property
    def json_serial(self):
        return self._json_serial

class RailRnaElasticAllJson(object):
    """ Constructs JSON for elastic mode + preprocess+align job flow. """
    def __init__(self, manifest, output_dir, intermediate_dir='./intermediate',
        force=False, aws_exe=None, profile='default', region='us-east-1',
        verbose=False, nucleotides_per_input=8000000, gzip_input=True,
        bowtie1_exe=None, bowtie1_idx='genome', bowtie1_build_exe=None,
        bowtie2_exe=None, bowtie2_build_exe=None, bowtie2_idx='genome',
        bowtie2_args='', samtools_exe=None, bedgraphtobigwig_exe=None,
        genome_partition_length=5000, max_readlet_size=25,
        readlet_config_size=32, min_readlet_size=15, readlet_interval=4,
        cap_size_multiplier=1.2, max_intron_size=500000, min_intron_size=10,
        min_exon_size=9, search_filter='none',
        motif_search_window_size=1000, max_gaps_mismatches=3,
        motif_radius=5, genome_bowtie1_args='-v 0 -a -m 80',
        transcriptome_bowtie2_args='-k 30', tie_margin=6, count_multiplier=6,
        normalize_percentile=0.75, very_replicable=False,
        do_not_output_bam_by_chr=False,
        output_sam=False, bam_basename='alignments', bed_basename='',
        log_uri=None, ami_version='3.2.1',
        visible_to_all_users=False, tags='',
        name='Rail-RNA Job Flow',
        action_on_failure='TERMINATE_JOB_FLOW',
        hadoop_jar=None,
        master_instance_count=1, master_instance_type='c1.xlarge',
        master_instance_bid_price=None, core_instance_count=1,
        core_instance_type=None, core_instance_bid_price=None,
        task_instance_count=0, task_instance_type=None,
        task_instance_bid_price=None, ec2_key_name=None, keep_alive=False,
        termination_protected=False, check_manifest=True,
        no_consistent_view=False, intermediate_lifetime=4):
        base = RailRnaErrors(manifest, output_dir, 
            intermediate_dir=intermediate_dir,
            force=force, aws_exe=aws_exe, profile=profile,
            region=region, verbose=verbose)
        RailRnaElastic(base, check_manifest=check_manifest, 
            log_uri=log_uri, ami_version=ami_version,
            visible_to_all_users=visible_to_all_users, tags=tags,
            name=name,
            action_on_failure=action_on_failure,
            hadoop_jar=hadoop_jar, master_instance_count=master_instance_count,
            master_instance_type=master_instance_type,
            master_instance_bid_price=master_instance_bid_price,
            core_instance_count=core_instance_count,
            core_instance_type=core_instance_type,
            core_instance_bid_price=core_instance_bid_price,
            task_instance_count=task_instance_count,
            task_instance_type=task_instance_type,
            task_instance_bid_price=task_instance_bid_price,
            ec2_key_name=ec2_key_name, keep_alive=keep_alive,
            termination_protected=termination_protected,
            no_consistent_view=no_consistent_view,
            intermediate_lifetime=intermediate_lifetime)
        RailRnaPreprocess(base, nucleotides_per_input=nucleotides_per_input,
            gzip_input=gzip_input)
        RailRnaAlign(base, elastic=True, bowtie1_exe=bowtie1_exe,
            bowtie1_idx=bowtie1_idx, bowtie1_build_exe=bowtie1_build_exe,
            bowtie2_exe=bowtie2_exe, bowtie2_build_exe=bowtie2_build_exe,
            bowtie2_idx=bowtie2_idx, bowtie2_args=bowtie2_args,
            samtools_exe=samtools_exe,
            bedgraphtobigwig_exe=bedgraphtobigwig_exe,
            genome_partition_length=genome_partition_length,
            max_readlet_size=max_readlet_size,
            readlet_config_size=readlet_config_size,
            min_readlet_size=min_readlet_size,
            readlet_interval=readlet_interval,
            cap_size_multiplier=cap_size_multiplier,
            max_intron_size=max_intron_size,
            min_intron_size=min_intron_size,
            search_filter=search_filter,
            min_exon_size=min_exon_size,
            motif_search_window_size=motif_search_window_size,
            max_gaps_mismatches=max_gaps_mismatches,
            motif_radius=motif_radius,
            genome_bowtie1_args=genome_bowtie1_args,
            transcriptome_bowtie2_args=transcriptome_bowtie2_args,
            count_multiplier=count_multiplier,
            tie_margin=tie_margin,
            normalize_percentile=normalize_percentile,
            very_replicable=very_replicable,
            do_not_output_bam_by_chr=do_not_output_bam_by_chr,
            output_sam=output_sam, bam_basename=bam_basename,
            bed_basename=bed_basename,
            s3_ansible=ab.S3Ansible(aws_exe=base.aws_exe,
                                        profile=base.profile))
        if base.errors:
            raise RuntimeError(
                    '\n'.join(
                            ['%d) %s' % (i+1, error) for i, error
                                in enumerate(base.errors)]
                        ) if len(base.errors) > 1 else base.errors[0]
                )
        self._json_serial = {}
        if base.core_instance_count > 0:
            reducer_count = base.core_instance_count \
                * base.instance_core_counts[base.core_instance_type]
        else:
            reducer_count = base.master_instance_count \
                * base.instance_core_counts[base.core_instance_type]
        prep_dir = path_join(True, base.intermediate_dir,
                                        'preprocess')
        push_dir = path_join(True, base.intermediate_dir,
                                        'preprocess', 'push')
        self._json_serial['Steps'] \
            = RailRnaElastic.hadoop_debugging_steps(base) + \
                steps(
                    RailRnaPreprocess.protosteps(base, prep_dir, push_dir,
                                                    elastic=True),
                    base.action_on_failure,
                    base.hadoop_jar, '/mnt/src/rna/steps',
                    reducer_count, base.intermediate_dir, unix=True,
                    no_consistent_view=base.no_consistent_view
                ) + \
                steps(
                    RailRnaAlign.protosteps(base, push_dir, elastic=True),
                    base.action_on_failure,
                    base.hadoop_jar, '/mnt/src/rna/steps',
                    reducer_count, base.intermediate_dir, unix=True,
                    no_consistent_view=base.no_consistent_view
                )
        self._json_serial['AmiVersion'] = base.ami_version
        if base.log_uri is not None:
            self._json_serial['LogUri'] = base.log_uri
        else:
            self._json_serial['LogUri'] = base.output_dir + '.logs'
        self._json_serial['Name'] = base.name
        self._json_serial['NewSupportedProducts'] = []
        self._json_serial['Tags'] = base.tags
        self._json_serial['VisibleToAllUsers'] = (
                'true' if base.visible_to_all_users else 'false'
            )
        self._json_serial['Instances'] = RailRnaElastic.instances(base)
        self._json_serial['BootstrapActions'] \
            = RailRnaAlign.bootstrap(base) \
            + RailRnaElastic.bootstrap(base)
        self.base = base
    
    @property
    def json_serial(self):
        return self._json_serial
