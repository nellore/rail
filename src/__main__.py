#!/usr/bin/env python
"""
rail-rna
Part of Rail-RNA

Main executable for Rail-RNA. Prints introductory message and "controls"
mode arguments {local, cloud} and flow arguments {preprocess, align, all}.
Command-line interface is inspired by git's and bowtie's.
"""

import sys
# Check for Python 2.7 immediately
if not (2 == sys.version_info[0] and sys.version_info[1] >= 7):
    print >>sys.stderr, ('Rail-RNA requires a version of Python >= 2.7 but '
                         '< 3.0. If an appropriate version of Python is '
                         'available at some path P that is not in PATH, try '
                         'running "P {0}". Otherwise, install an '
                         'appropriate version of Python and rerun.').format(
                                ' '.join(sys.argv)
                            )
    sys.exit(1)
import os
base_path = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
driver_path = os.path.join(base_path, 'rna', 'driver')
import site
site.addsitedir(driver_path)

'''Is this the installer? Assumed to be if Rail-RNA's zipped up. Check if the
currently executing script is in a zip.'''
import zipfile
containing_dir = os.path.abspath(os.path.dirname(__file__))
if zipfile.is_zipfile(containing_dir):
    # In a zip; use install parser
    import argparse
    import exe_paths
    from rna_config import rail_help_wrapper
    parser = argparse.ArgumentParser(
                    formatter_class=rail_help_wrapper,
                    add_help=True
                )
    parser.add_argument('-i', '--install-dir', type=str, required=False,
            metavar='<dir>',
            default=None,
            help=('directory in which to install Rail-RNA (def: '
                  '"/usr/local/raildotbio" if installing for all users and '
                  '"~/raildotbio" if installing for current user)')
        )
    parser.add_argument('-n', '--no-dependencies', action='store_const',
            const=True,
            default=False,
            help='installs Rail-RNA without any of its dependencies'
        )
    parser.add_argument('-p', '--prep-dependencies', action='store_const',
            const=True,
            default=False,
            help=('installs Rail-RNA with only dependencies required for its '
                  'preprocess job flow; overrided by --no-dependencies')
        )
    parser.add_argument('-y', '--yes', action='store_const',
            const=True,
            default=False,
            help=('answers "yes" to all user prompts, installing for all '
                  'users; overrided by --me')
        )
    parser.add_argument('-m', '--me', action='store_const',
            const=True,
            default=False,
            help=('answers "no" to the user prompt "Install for all users?" '
                  'and "yes" to all other user prompts')
        )
    parser.add_argument('-s', '--symlink-dependencies', action='store_const',
            const=True,
            default=False,
            help=('symlinks all installed dependencies to /usr/local/bin if '
                  'installing for all users')
        )
    parser.add_argument('--curl', type=str, required=False, metavar='<exe>',
            default=exe_paths.curl,
            help=('path to cURL executable (def: %s)'
                    % (exe_paths.curl if exe_paths.curl is not None
                        else 'curl')))
    args = parser.parse_args()
    from rna_installer import RailRnaInstaller
    with RailRnaInstaller(containing_dir, curl_exe=args.curl,
                            install_dir=args.install_dir,
                            no_dependencies=args.no_dependencies,
                            prep_dependencies=args.prep_dependencies,
                            add_symlinks=args.symlink_dependencies,
                            yes=args.yes, me=args.me) as railrna_installer:
        railrna_installer.install()
    sys.exit(0)

site.addsitedir(base_path)
from rna_config import *
from rna_config import _warning_message, _executable
from dooplicity.tools import which
import json
import subprocess
from argparse import SUPPRESS
from version import version_number
import datetime

_usage_message = \
u"""rail-rna <job flow> <mode> <[args]>

  <job flow>       {{prep, align, go}}
                     prep: preprocess reads listed in a required manifest
                       file (specified with --manifest)
                     align: align preprocessed reads (specified with --input)
                     go: perform prep and align in succession
  <mode>           {{local, parallel, elastic}}
                     local: run Rail-RNA on this computer
                     parallel: run Rail-RNA on all active IPython engines
                     elastic: run Rail-RNA on Amazon Elastic MapReduce.
                       Requires that the user sign up for Amazon Web Services

{0} Rail-RNA v{1} by Abhi Nellore (anellore@jhu.edu; nellore.github.io)

Rail-RNA is a scalable MapReduce pipeline that can analyze many RNA-seq
samples at once. To view help for a given combination of <job flow> and
<mode>, specify both, then add -h/--help.""".format(u'\u2200', version_number)

if sys.stdout.encoding != 'UTF-8':
    # Correct for consoles that don't support UTF-8
    import unicodedata
    _usage_message = unicodedata.normalize(
                            'NFKD', _usage_message
                        ).encode('ascii', 'ignore')

class Launcher(object):
    """ Facilitates replacing the current process with a Dooplicity runner. """

    def __init__(self, force=False, num_processes=1, keep_intermediates=False,
                    gzip_intermediates=False, gzip_level=3,
                    sort_memory_cap=0.2, max_task_attempts=4,
                    region='us-east-1', log=None, scratch=None,
                    ipython_profile=None, ipcontroller_json=None, common=None,
                    direct_write=False, json=False, sort=None,
                    profile=None):
        self.force = force
        self.num_processes = num_processes
        self.keep_intermediates = keep_intermediates
        self.gzip_intermediates = gzip_intermediates
        self.gzip_level = gzip_level
        self.sort_memory_cap = sort_memory_cap
        self.max_task_attempts = max_task_attempts
        self.region = region
        self.log = log
        self.scratch = scratch
        self.direct_write = direct_write
        self.ipython_profile = ipython_profile
        self.ipcontroller_json = ipcontroller_json
        self.common = common
        self.json = json
        self.sort = sort
        self.profile = profile

    def run(self, mode, payload):
        """ Replaces current process, using PyPy if it's available.

            Transferring control allows Dooplicity to handle errors from here
            on out.

            mode: 'local' or 'elastic'; launches emr_simulator.py or 
                emr_runner.py, respectively
            payload: string with json payload to copy to stdin
                of replacement process
        """
        if self.json:
            print json.dumps(json.loads(payload), sort_keys=True,
                             indent=4, separators=(',', ': '))
            quit()
        #temp=json.loads(payload)
        #del temp['ServiceRole']
        #del temp['JobFlowRole']
        #payload = json.dumps(temp)
        read_pipe, write_pipe = os.pipe()
        if os.fork() != 0:
            # Parent process; read from child after determining executable
            if mode == 'local':
                print_to_screen(_warning_message)
                runner_args = [_executable, os.path.join(
                                                    base_path,
                                                    'dooplicity',
                                                    'emr_simulator.py'
                                                ),
                                '-p', str(self.num_processes),
                                '-b', os.path.join(base_path, 
                                        'rna', 'driver', 'rail-rna.txt'),
                                '--memcap', str(self.sort_memory_cap),
                                '--max-attempts', str(self.max_task_attempts)]
                if self.sort:
                    runner_args.extend(['--sort', self.sort])
                if self.force:
                    runner_args.append('-f')
                if self.keep_intermediates:
                    runner_args.append('--keep-intermediates')
                if self.gzip_intermediates:
                    runner_args.extend(['--gzip-outputs', '--gzip-level',
                                            str(self.gzip_level)])
                if self.log:
                    runner_args.extend(['-l', os.path.abspath(self.log)])
                os.dup2(read_pipe, sys.stdin.fileno())
                os.close(read_pipe)
                os.close(write_pipe)
                os.execv(_executable, runner_args)
            elif mode == 'parallel':
                print_to_screen('Launching Dooplicity runner with Python...')
                parallel_executable = which('python')
                # sys.executable had better find IPython
                runner_args = [parallel_executable, os.path.join(
                                                    base_path,
                                                    'dooplicity',
                                                    'emr_simulator.py'
                                                ),
                                '-b', os.path.join(base_path, 
                                        'rna', 'driver', 'rail-rna.txt'),
                                '--ipy',
                                '--memcap', str(self.sort_memory_cap),
                                '--max-attempts', str(self.max_task_attempts)]
                if self.sort:
                    runner_args.extend(['--sort', self.sort])
                if self.force:
                    runner_args.append('-f')
                if self.keep_intermediates:
                    runner_args.append('--keep-intermediates')
                if self.direct_write:
                    runner_args.append('--direct-write')
                if self.gzip_intermediates:
                    runner_args.extend(['--gzip-outputs', '--gzip-level',
                                            str(self.gzip_level)])
                if self.log:
                    runner_args.extend(['-l', os.path.abspath(self.log)])
                if self.common:
                    runner_args.extend(['--common',
                        os.path.abspath(self.common)])
                if self.scratch:
                    runner_args.extend(['--scratch', self.scratch])
                else:
                    runner_args.extend(['--scratch', '-'])
                if self.ipython_profile:
                    runner_args.extend(['--ipy-profile', self.ipython_profile])
                if self.ipcontroller_json:
                    runner_args.extend(['--ipcontroller-json',
                                            self.ipcontroller_json])
                os.dup2(read_pipe, sys.stdin.fileno())
                os.close(read_pipe)
                os.close(write_pipe)
                os.execv(sys.executable, runner_args)
            else:
                runner_args = [_executable, os.path.join(
                                                    base_path,
                                                    'dooplicity',
                                                    'emr_runner.py'
                                                ),
                                '-b', os.path.join(base_path, 
                                        'rna', 'driver', 'rail-rna.txt')]
                if self.force:
                    runner_args.append('-f')
                if self.profile:
                    runner_args.extend(['--profile', self.profile])
                runner_args.extend(['-r', self.region])
                os.dup2(read_pipe, sys.stdin.fileno())
                os.close(read_pipe)
                os.close(write_pipe)
                os.execv(_executable, runner_args)
        else:
            os.write(write_pipe, payload)
            os.close(write_pipe)
            ###SCRIPT TERMINATES HERE###

if __name__ == '__main__':
    parser = RailParser(
            usage=_usage_message,
            add_help=False
        )
    flow_parsers = parser.add_subparsers(
            dest='job_flow'
        )
    prep_mode_parser = flow_parsers.add_parser('prep', add_help=False,
                                                    usage=_usage_message)
    align_mode_parser = flow_parsers.add_parser('align', add_help=False,
                                                    usage=_usage_message)
    go_mode_parser = flow_parsers.add_parser('go', add_help=False,
                                                    usage=_usage_message)
    prep_parsers = prep_mode_parser.add_subparsers(dest='prep_mode')
    prep_local_parser = prep_parsers.add_parser(
                                    'local',
                                    usage=general_usage('prep local',
                                        '-m <file> '),
                                    formatter_class=rail_help_wrapper,
                                    add_help=False
                                )
    prep_parallel_parser = prep_parsers.add_parser(
                                    'parallel',
                                    usage=general_usage('prep parallel',
                                        '-m <file> '),
                                    formatter_class=rail_help_wrapper,
                                    add_help=False
                                )
    prep_elastic_parser = prep_parsers.add_parser(
                                    'elastic',
                                    usage=general_usage('prep elastic',
                                        '-m <file> -o <s3_dir> -c <int> '),
                                    formatter_class=rail_help_wrapper,
                                    add_help=False
                                )
    align_parsers = align_mode_parser.add_subparsers(dest='align_mode')
    align_local_parser = align_parsers.add_parser(
                                    'local',
                                    usage=general_usage('align local',
                                        '-m <file> -i <dir> '
                                        '-x <idx | idx,idx> '),
                                    formatter_class=rail_help_wrapper,
                                    add_help=False
                                )
    align_parallel_parser = align_parsers.add_parser(
                                    'parallel',
                                    usage=general_usage('align parallel',
                                        '-m <file> -i <dir> '
                                        '-x <idx | idx,idx> '),
                                    formatter_class=rail_help_wrapper,
                                    add_help=False
                                )
    align_elastic_parser = align_parsers.add_parser(
                                    'elastic',
                                    usage=general_usage('align elastic',
                                        '-m <file> -i <s3_dir> -a '
                                        '<choice | tgz> \r\n       '
                                        '-c <int> -o <s3_dir> '),
                                    formatter_class=rail_help_wrapper,
                                    add_help=False
                                )
    go_parsers = go_mode_parser.add_subparsers(dest='go_mode')
    go_local_parser = go_parsers.add_parser(
                                    'local',
                                    usage=general_usage('go local',
                                        '-m <file> -x <idx | idx,idx> '),
                                    formatter_class=rail_help_wrapper,
                                    add_help=False
                                )
    go_parallel_parser = go_parsers.add_parser(
                                    'parallel',
                                    usage=general_usage('go parallel',
                                        '-m <file> -x <idx | idx,idx> '),
                                    formatter_class=rail_help_wrapper,
                                    add_help=False
                                )
    go_elastic_parser = go_parsers.add_parser(
                                    'elastic',
                                    usage=general_usage('go elastic',
                                        '-m <file> -o <s3_dir> '
                                        '-a <choice/tgz> \r\n       '
                                        '-c <int> '),
                                    formatter_class=rail_help_wrapper,
                                    add_help=False
                                )
    '''Group arguments; code seems more involved than it is. Could perhaps be
    condensed, but each mode+job flow could ultimately become more
    customized, rendering more general code unusable.'''
    prep_local_required \
        = prep_local_parser.add_argument_group('required arguments')
    prep_local_output \
        = prep_local_parser.add_argument_group('output options')
    prep_local_general \
        = prep_local_parser.add_argument_group('general options')
    prep_local_exec \
        = prep_local_parser.add_argument_group('dependencies')
    prep_parallel_required \
        = prep_parallel_parser.add_argument_group('required arguments')
    prep_parallel_output \
        = prep_parallel_parser.add_argument_group('output options')
    prep_parallel_general \
        = prep_parallel_parser.add_argument_group('general options')
    prep_parallel_exec \
        = prep_parallel_parser.add_argument_group('dependencies')
    prep_elastic_required \
        = prep_elastic_parser.add_argument_group('required arguments')
    prep_elastic_output \
        = prep_elastic_parser.add_argument_group('output options')
    prep_elastic_general \
        = prep_elastic_parser.add_argument_group('general options')
    prep_elastic_exec \
        = prep_elastic_parser.add_argument_group('dependencies')
    prep_elastic_details \
        = prep_elastic_parser.add_argument_group('Elastic MapReduce options')
    align_local_required \
        = align_local_parser.add_argument_group('required arguments')
    align_local_output \
        = align_local_parser.add_argument_group('output options')
    align_local_general \
        = align_local_parser.add_argument_group('general options')
    align_local_exec \
        = align_local_parser.add_argument_group('dependencies')
    align_local_algo \
        = align_local_parser.add_argument_group('algorithm options')
    align_parallel_required \
        = align_parallel_parser.add_argument_group('required arguments')
    align_parallel_output \
        = align_parallel_parser.add_argument_group('output options')
    align_parallel_general \
        = align_parallel_parser.add_argument_group('general options')
    align_parallel_exec \
        = align_parallel_parser.add_argument_group('dependencies')
    align_parallel_algo \
        = align_parallel_parser.add_argument_group('algorithm options')
    align_elastic_required \
        = align_elastic_parser.add_argument_group('required arguments')
    align_elastic_output \
        = align_elastic_parser.add_argument_group('output options')
    align_elastic_general \
        = align_elastic_parser.add_argument_group('general options')
    align_elastic_exec \
        = align_elastic_parser.add_argument_group('dependencies')
    align_elastic_details \
        = align_elastic_parser.add_argument_group('Elastic MapReduce options')
    align_elastic_algo \
        = align_elastic_parser.add_argument_group('algorithm options')
    go_local_required \
        = go_local_parser.add_argument_group('required arguments')
    go_local_output \
        = go_local_parser.add_argument_group('output options')
    go_local_general \
        = go_local_parser.add_argument_group('general options')
    go_local_exec \
        = go_local_parser.add_argument_group('dependencies')
    go_local_algo \
        = go_local_parser.add_argument_group('algorithm options')
    go_parallel_required \
        = go_parallel_parser.add_argument_group('required arguments')
    go_parallel_output \
        = go_parallel_parser.add_argument_group('output options')
    go_parallel_general \
        = go_parallel_parser.add_argument_group('general options')
    go_parallel_exec \
        = go_parallel_parser.add_argument_group('dependencies')
    go_parallel_algo \
        = go_parallel_parser.add_argument_group('algorithm options')
    go_elastic_required \
        = go_elastic_parser.add_argument_group('required arguments')
    go_elastic_output \
        = go_elastic_parser.add_argument_group('output options')
    go_elastic_general \
        = go_elastic_parser.add_argument_group('general options')
    go_elastic_exec \
        = go_elastic_parser.add_argument_group('dependencies')
    go_elastic_details \
        = go_elastic_parser.add_argument_group('Elastic MapReduce options')
    go_elastic_algo \
        = go_elastic_parser.add_argument_group('algorithm options')
    # Add helps manually to general options
    for subparser in [prep_local_general, align_local_general,
                        go_local_general, prep_parallel_general,
                        align_parallel_general, go_parallel_general,
                        prep_elastic_general, align_elastic_general,
                        go_elastic_general]:
        subparser.add_argument(
                    '-h', '--help',
                    action='help', default=SUPPRESS,
                    help='show this help message and exit'
                )
        subparser.add_argument(
                    '-v', '--version',
                    action='version',
                    version=('Rail-RNA v{0}'.format(version_number)),
                    help='show version information and exit'
                )
        subparser.add_argument(
                    '-j', '--json',
                    action='store_const', const=True,
                    default=False,
                    help=('print job flow JSON to stdout and exit')
                )
    parser.add_argument(
            '-v', '--version',
            action='version',
            version=('Rail-RNA v{0}'.format(version_number))
        )
    RailRnaErrors.add_args(general_parser=go_elastic_general,
                            exec_parser=go_elastic_exec,
                            required_parser=go_elastic_required)
    RailRnaErrors.add_args(general_parser=go_local_general,
                            exec_parser=go_local_exec,
                            required_parser=go_local_required)
    RailRnaErrors.add_args(general_parser=go_parallel_general,
                            exec_parser=go_parallel_exec,
                            required_parser=go_parallel_required)
    RailRnaErrors.add_args(general_parser=align_elastic_general,
                            exec_parser=align_elastic_exec,
                            required_parser=align_elastic_required)
    RailRnaErrors.add_args(general_parser=align_local_general,
                            exec_parser=align_local_exec,
                            required_parser=align_local_required)
    RailRnaErrors.add_args(general_parser=align_parallel_general,
                            exec_parser=align_parallel_exec,
                            required_parser=align_parallel_required)
    RailRnaErrors.add_args(general_parser=prep_elastic_general,
                            exec_parser=prep_elastic_exec,
                            required_parser=prep_elastic_required)
    RailRnaErrors.add_args(general_parser=prep_local_general,
                            exec_parser=prep_local_exec,
                            required_parser=prep_local_required)
    RailRnaLocal.add_args(required_parser=go_local_required,
                            general_parser=go_local_general,
                            output_parser=go_local_output,
                            exec_parser=go_local_exec,
                            prep=False, align=False)
    RailRnaLocal.add_args(required_parser=align_local_required,
                            general_parser=align_local_general,
                            output_parser=align_local_output,
                            exec_parser=align_local_exec,
                            prep=False, align=True)
    RailRnaLocal.add_args(required_parser=prep_local_required,
                            general_parser=prep_local_general,
                            output_parser=prep_local_output,
                            exec_parser=prep_local_exec,
                            prep=True, align=False)
    RailRnaErrors.add_args(general_parser=prep_parallel_general,
                            exec_parser=prep_parallel_exec,
                            required_parser=prep_parallel_required)
    RailRnaLocal.add_args(required_parser=go_parallel_required,
                            general_parser=go_parallel_general,
                            output_parser=go_parallel_output,
                            exec_parser=go_parallel_exec,
                            prep=False, align=False, parallel=True)
    RailRnaLocal.add_args(required_parser=align_parallel_required,
                            general_parser=align_parallel_general,
                            output_parser=align_parallel_output,
                            exec_parser=align_parallel_exec,
                            prep=False, align=True, parallel=True)
    RailRnaLocal.add_args(required_parser=prep_parallel_required,
                            general_parser=prep_parallel_general,
                            output_parser=prep_parallel_output,
                            exec_parser=prep_parallel_exec,
                            prep=True, align=False, parallel=True)
    RailRnaElastic.add_args(required_parser=go_elastic_required,
                            general_parser=go_elastic_general,
                            output_parser=go_elastic_output,
                            elastic_parser=go_elastic_details,
                            align=False)
    RailRnaElastic.add_args(required_parser=align_elastic_required,
                            general_parser=align_elastic_general,
                            output_parser=align_elastic_output,
                            elastic_parser=align_elastic_details,
                            align=True)
    RailRnaElastic.add_args(required_parser=prep_elastic_required,
                            general_parser=prep_elastic_general,
                            output_parser=prep_elastic_output,
                            elastic_parser=prep_elastic_details,
                            align=False)
    RailRnaPreprocess.add_args(general_parser=prep_elastic_general,
                               output_parser=prep_elastic_output,
                               elastic=True)
    RailRnaPreprocess.add_args(general_parser=prep_local_general,
                               output_parser=prep_local_output,
                               elastic=False)
    RailRnaPreprocess.add_args(general_parser=prep_parallel_general,
                               output_parser=prep_parallel_output,
                               elastic=False)
    RailRnaPreprocess.add_args(general_parser=go_elastic_general,
                               output_parser=go_elastic_general,
                               elastic=True)
    RailRnaPreprocess.add_args(general_parser=go_local_general,
                               output_parser=go_local_general,
                               elastic=False)
    RailRnaPreprocess.add_args(general_parser=go_parallel_general,
                               output_parser=go_parallel_general,
                               elastic=False)
    RailRnaAlign.add_args(required_parser=align_local_required,
                          exec_parser=align_local_exec,
                          output_parser=align_local_output,
                          algo_parser=align_local_algo, elastic=False)
    RailRnaAlign.add_args(required_parser=align_parallel_required,
                          exec_parser=align_parallel_exec,
                          output_parser=align_parallel_output,
                          algo_parser=align_parallel_algo, elastic=False)
    RailRnaAlign.add_args(required_parser=align_elastic_required,
                          exec_parser=align_elastic_exec,
                          output_parser=align_elastic_output,
                          algo_parser=align_elastic_algo, elastic=True)
    RailRnaAlign.add_args(required_parser=go_local_required,
                          exec_parser=go_local_exec,
                          output_parser=go_local_output,
                          algo_parser=go_local_algo, elastic=False)
    RailRnaAlign.add_args(required_parser=go_parallel_required,
                          exec_parser=go_parallel_exec,
                          output_parser=go_parallel_output,
                          algo_parser=go_parallel_algo, elastic=False)
    RailRnaAlign.add_args(required_parser=go_elastic_required,
                          exec_parser=go_elastic_exec,
                          output_parser=go_elastic_output,
                          algo_parser=go_elastic_algo, elastic=True)
    args = parser.parse_args()
    print_to_screen('Loading...', newline=True, carriage_return=False)
    if args.job_flow == 'go' and args.go_mode == 'local':
        mode = 'local'
        json_creator = RailRnaLocalAllJson(
                args.manifest, args.output,
                isofrag_idx=args.isofrag_idx,
                intermediate_dir=args.log,
                force=args.force, aws_exe=args.aws, profile=args.profile,
                verbose=args.verbose,
                nucleotides_per_input=args.nucleotides_per_input,
                gzip_input=(not args.do_not_gzip_input),
                do_not_bin_quals=args.do_not_bin_quals,
                short_read_names=args.short_read_names,
                skip_bad_records=args.skip_bad_records,
                ignore_missing_sra_samples=args.ignore_missing_sra_samples,
                bowtie_idx=args.bowtie_idx,
                bowtie1_exe=args.bowtie1, bowtie2_exe=args.bowtie2,
                bowtie1_build_exe=args.bowtie1_build,
                bowtie2_build_exe=args.bowtie2_build,
                k=args.k, bowtie2_args=args.bowtie2_args,
                samtools_exe=args.samtools,
                bedgraphtobigwig_exe=args.bedgraphtobigwig,
                partition_length=args.partition_length,
                max_readlet_size=args.max_readlet_size,
                readlet_config_size=args.readlet_config_size,
                min_readlet_size=args.min_readlet_size,
                readlet_interval=args.readlet_interval,
                cap_size_multiplier=args.cap_size_multiplier,
                max_intron_size=args.max_intron_size,
                min_intron_size=args.min_intron_size,
                min_exon_size=args.min_exon_size,
                library_size=args.library_size,
                search_filter=args.search_filter,
                motif_search_window_size=args.motif_search_window_size,
                max_gaps_mismatches=args.max_gaps_mismatches,
                motif_radius=args.motif_radius,
                genome_bowtie1_args=args.genome_bowtie1_args,
                junction_criteria=args.junction_criteria,
                indel_criteria=args.indel_criteria,
                transcriptome_bowtie2_args=args.transcriptome_bowtie2_args,
                experimental=args.experimental,
                count_multiplier=args.count_multiplier,
                max_refs_per_strand=args.max_refs_per_strand,
                tie_margin=args.tie_margin,
                normalize_percentile=args.normalize_percentile,
                transcriptome_indexes_per_sample=\
                    args.transcriptome_indexes_per_sample,
                drop_deletions=args.drop_deletions,
                do_not_output_bam_by_chr=args.do_not_output_bam_by_chr,
                do_not_output_ave_bw_by_chr=args.do_not_output_ave_bw_by_chr,
                do_not_drop_polyA_tails=args.do_not_drop_polyA_tails,
                deliverables=args.deliverables,
                bam_basename=args.bam_basename,
                tsv_basename=args.tsv_basename,
                bed_basename=args.bed_basename,
                num_processes=args.num_processes,
                gzip_intermediates=args.gzip_intermediates,
                gzip_level=args.gzip_level,
                sort_memory_cap=args.sort_memory_cap,
                max_task_attempts=args.max_task_attempts,
                keep_intermediates=args.keep_intermediates,
                check_manifest=(not args.do_not_check_manifest),
                sort_exe=args.sort,
                scratch=args.scratch,
                fastq_dump_exe=args.fastq_dump,
                vdb_config_exe=args.vdb_config
            )
    elif args.job_flow == 'align' and args.align_mode == 'local':
        mode = 'local'
        json_creator = RailRnaLocalAlignJson(
                args.manifest, args.output, args.input,
                isofrag_idx=args.isofrag_idx,
                intermediate_dir=args.log,
                force=args.force, aws_exe=args.aws, profile=args.profile,
                verbose=args.verbose,
                bowtie_idx=args.bowtie_idx,
                bowtie1_exe=args.bowtie1, bowtie2_exe=args.bowtie2,
                bowtie1_build_exe=args.bowtie1_build,
                bowtie2_build_exe=args.bowtie2_build,
                k=args.k, bowtie2_args=args.bowtie2_args,
                samtools_exe=args.samtools,
                bedgraphtobigwig_exe=args.bedgraphtobigwig,
                partition_length=args.partition_length,
                max_readlet_size=args.max_readlet_size,
                readlet_config_size=args.readlet_config_size,
                min_readlet_size=args.min_readlet_size,
                readlet_interval=args.readlet_interval,
                cap_size_multiplier=args.cap_size_multiplier,
                max_intron_size=args.max_intron_size,
                min_intron_size=args.min_intron_size,
                min_exon_size=args.min_exon_size,
                library_size=args.library_size,
                search_filter=args.search_filter,
                motif_search_window_size=args.motif_search_window_size,
                max_gaps_mismatches=args.max_gaps_mismatches,
                motif_radius=args.motif_radius,
                genome_bowtie1_args=args.genome_bowtie1_args,
                junction_criteria=args.junction_criteria,
                indel_criteria=args.indel_criteria,
                transcriptome_bowtie2_args=args.transcriptome_bowtie2_args,
                experimental=args.experimental,
                count_multiplier=args.count_multiplier,
                max_refs_per_strand=args.max_refs_per_strand,
                tie_margin=args.tie_margin,
                normalize_percentile=args.normalize_percentile,
                transcriptome_indexes_per_sample=\
                    args.transcriptome_indexes_per_sample,
                drop_deletions=args.drop_deletions,
                do_not_output_bam_by_chr=args.do_not_output_bam_by_chr,
                do_not_output_ave_bw_by_chr=args.do_not_output_ave_bw_by_chr,
                do_not_drop_polyA_tails=args.do_not_drop_polyA_tails,
                deliverables=args.deliverables,
                bam_basename=args.bam_basename,
                tsv_basename=args.tsv_basename,
                bed_basename=args.bed_basename,
                num_processes=args.num_processes,
                gzip_intermediates=args.gzip_intermediates,
                gzip_level=args.gzip_level,
                sort_memory_cap=args.sort_memory_cap,
                max_task_attempts=args.max_task_attempts,
                keep_intermediates=args.keep_intermediates,
                sort_exe=args.sort,
                scratch=args.scratch
            )
    elif args.job_flow == 'prep' and args.prep_mode == 'local':
        mode = 'local'
        json_creator = RailRnaLocalPreprocessJson(
                args.manifest, args.output,
                intermediate_dir=args.log,
                force=args.force, aws_exe=args.aws, profile=args.profile,
                verbose=args.verbose,
                nucleotides_per_input=args.nucleotides_per_input,
                gzip_input=(not args.do_not_gzip_input),
                do_not_bin_quals=args.do_not_bin_quals,
                short_read_names=args.short_read_names,
                skip_bad_records=args.skip_bad_records,
                ignore_missing_sra_samples=args.ignore_missing_sra_samples,
                num_processes=args.num_processes,
                gzip_intermediates=args.gzip_intermediates,
                gzip_level=args.gzip_level,
                sort_memory_cap=args.sort_memory_cap,
                max_task_attempts=args.max_task_attempts,
                keep_intermediates=args.keep_intermediates,
                check_manifest=(not args.do_not_check_manifest),
                sort_exe=args.sort,
                scratch=args.scratch,
                fastq_dump_exe=args.fastq_dump,
                vdb_config_exe=args.vdb_config
            )
    elif args.job_flow == 'go' and args.go_mode == 'parallel':
        mode = 'parallel'
        json_creator = RailRnaParallelAllJson(
                args.manifest, args.output,
                isofrag_idx=args.isofrag_idx,
                intermediate_dir=args.log,
                force=args.force, aws_exe=args.aws, profile=args.profile,
                verbose=args.verbose,
                nucleotides_per_input=args.nucleotides_per_input,
                gzip_input=(not args.do_not_gzip_input),
                do_not_bin_quals=args.do_not_bin_quals,
                short_read_names=args.short_read_names,
                skip_bad_records=args.skip_bad_records,
                ignore_missing_sra_samples=args.ignore_missing_sra_samples,
                bowtie_idx=args.bowtie_idx,
                bowtie1_exe=args.bowtie1, bowtie2_exe=args.bowtie2,
                bowtie1_build_exe=args.bowtie1_build,
                bowtie2_build_exe=args.bowtie2_build,
                k=args.k, bowtie2_args=args.bowtie2_args,
                samtools_exe=args.samtools,
                bedgraphtobigwig_exe=args.bedgraphtobigwig,
                partition_length=args.partition_length,
                max_readlet_size=args.max_readlet_size,
                readlet_config_size=args.readlet_config_size,
                min_readlet_size=args.min_readlet_size,
                readlet_interval=args.readlet_interval,
                cap_size_multiplier=args.cap_size_multiplier,
                max_intron_size=args.max_intron_size,
                min_intron_size=args.min_intron_size,
                min_exon_size=args.min_exon_size,
                library_size=args.library_size,
                search_filter=args.search_filter,
                motif_search_window_size=args.motif_search_window_size,
                max_gaps_mismatches=args.max_gaps_mismatches,
                motif_radius=args.motif_radius,
                genome_bowtie1_args=args.genome_bowtie1_args,
                junction_criteria=args.junction_criteria,
                indel_criteria=args.indel_criteria,
                transcriptome_bowtie2_args=args.transcriptome_bowtie2_args,
                experimental=args.experimental,
                count_multiplier=args.count_multiplier,
                max_refs_per_strand=args.max_refs_per_strand,
                tie_margin=args.tie_margin,
                normalize_percentile=args.normalize_percentile,
                transcriptome_indexes_per_sample=\
                    args.transcriptome_indexes_per_sample,
                drop_deletions=args.drop_deletions,
                do_not_output_bam_by_chr=args.do_not_output_bam_by_chr,
                do_not_output_ave_bw_by_chr=args.do_not_output_ave_bw_by_chr,
                do_not_drop_polyA_tails=args.do_not_drop_polyA_tails,
                deliverables=args.deliverables,
                bam_basename=args.bam_basename,
                tsv_basename=args.tsv_basename,
                bed_basename=args.bed_basename,
                gzip_intermediates=args.gzip_intermediates,
                gzip_level=args.gzip_level,
                sort_memory_cap=args.sort_memory_cap,
                max_task_attempts=args.max_task_attempts,
                keep_intermediates=args.keep_intermediates,
                check_manifest=(not args.do_not_check_manifest),
                ipython_profile=args.ipython_profile,
                ipcontroller_json=args.ipcontroller_json,
                scratch=args.scratch,
                direct_write=args.direct_write,
                do_not_copy_index_to_nodes=args.do_not_copy_index_to_nodes,
                sort_exe=args.sort,
                fastq_dump_exe=args.fastq_dump,
                vdb_config_exe=args.vdb_config
            )
    elif args.job_flow == 'align' and args.align_mode == 'parallel':
        mode = 'parallel'
        json_creator = RailRnaParallelAlignJson(
                args.manifest, args.output, args.input,
                isofrag_idx=args.isofrag_idx,
                intermediate_dir=args.log,
                force=args.force, aws_exe=args.aws, profile=args.profile,
                verbose=args.verbose,
                bowtie_idx=args.bowtie_idx,
                bowtie1_exe=args.bowtie1, bowtie2_exe=args.bowtie2,
                bowtie1_build_exe=args.bowtie1_build,
                bowtie2_build_exe=args.bowtie2_build,
                k=args.k, bowtie2_args=args.bowtie2_args,
                samtools_exe=args.samtools,
                bedgraphtobigwig_exe=args.bedgraphtobigwig,
                partition_length=args.partition_length,
                max_readlet_size=args.max_readlet_size,
                readlet_config_size=args.readlet_config_size,
                min_readlet_size=args.min_readlet_size,
                readlet_interval=args.readlet_interval,
                cap_size_multiplier=args.cap_size_multiplier,
                max_intron_size=args.max_intron_size,
                min_intron_size=args.min_intron_size,
                min_exon_size=args.min_exon_size,
                library_size=args.library_size,
                search_filter=args.search_filter,
                motif_search_window_size=args.motif_search_window_size,
                max_gaps_mismatches=args.max_gaps_mismatches,
                motif_radius=args.motif_radius,
                genome_bowtie1_args=args.genome_bowtie1_args,
                junction_criteria=args.junction_criteria,
                indel_criteria=args.indel_criteria,
                transcriptome_bowtie2_args=args.transcriptome_bowtie2_args,
                experimental=args.experimental,
                count_multiplier=args.count_multiplier,
                max_refs_per_strand=args.max_refs_per_strand,
                tie_margin=args.tie_margin,
                normalize_percentile=args.normalize_percentile,
                transcriptome_indexes_per_sample=\
                    args.transcriptome_indexes_per_sample,
                drop_deletions=args.drop_deletions,
                do_not_output_bam_by_chr=args.do_not_output_bam_by_chr,
                do_not_output_ave_bw_by_chr=args.do_not_output_ave_bw_by_chr,
                do_not_drop_polyA_tails=args.do_not_drop_polyA_tails,
                deliverables=args.deliverables,
                bam_basename=args.bam_basename,
                tsv_basename=args.tsv_basename,
                bed_basename=args.bed_basename,
                gzip_intermediates=args.gzip_intermediates,
                gzip_level=args.gzip_level,
                sort_memory_cap=args.sort_memory_cap,
                max_task_attempts=args.max_task_attempts,
                keep_intermediates=args.keep_intermediates,
                ipython_profile=args.ipython_profile,
                ipcontroller_json=args.ipcontroller_json,
                scratch=args.scratch,
                direct_write=args.direct_write,
                do_not_copy_index_to_nodes=args.do_not_copy_index_to_nodes,
                sort_exe=args.sort
            )
    elif args.job_flow == 'prep' and args.prep_mode == 'parallel':
        mode = 'parallel'
        json_creator = RailRnaParallelPreprocessJson(
                args.manifest, args.output,
                intermediate_dir=args.log,
                force=args.force, aws_exe=args.aws, profile=args.profile,
                verbose=args.verbose,
                nucleotides_per_input=args.nucleotides_per_input,
                gzip_input=(not args.do_not_gzip_input),
                do_not_bin_quals=args.do_not_bin_quals,
                short_read_names=args.short_read_names,
                skip_bad_records=args.skip_bad_records,
                ignore_missing_sra_samples=args.ignore_missing_sra_samples,
                gzip_intermediates=args.gzip_intermediates,
                gzip_level=args.gzip_level,
                sort_memory_cap=args.sort_memory_cap,
                max_task_attempts=args.max_task_attempts,
                keep_intermediates=args.keep_intermediates,
                check_manifest=(not args.do_not_check_manifest),
                ipython_profile=args.ipython_profile,
                ipcontroller_json=args.ipcontroller_json,
                scratch=args.scratch,
                direct_write=args.direct_write,
                sort_exe=args.sort,
                fastq_dump_exe=args.fastq_dump,
                vdb_config_exe=args.vdb_config
            )
    elif args.job_flow == 'go' and args.go_mode == 'elastic':
        mode = 'elastic'
        json_creator = RailRnaElasticAllJson(
                args.manifest, args.output,
                assembly=args.assembly,
                do_not_bin_quals=args.do_not_bin_quals,
                short_read_names=args.short_read_names,
                skip_bad_records=args.skip_bad_records,
                ignore_missing_sra_samples=args.ignore_missing_sra_samples,
                isofrag_idx=args.isofrag_idx,
                intermediate_dir=args.intermediate,
                force=args.force, aws_exe=args.aws, profile=args.profile,
                region=args.region,
                service_role=args.service_role,
                instance_profile=args.instance_profile,
                verbose=args.verbose,
                k=args.k, bowtie2_args=args.bowtie2_args,
                partition_length=args.partition_length,
                max_readlet_size=args.max_readlet_size,
                readlet_config_size=args.readlet_config_size,
                min_readlet_size=args.min_readlet_size,
                readlet_interval=args.readlet_interval,
                cap_size_multiplier=args.cap_size_multiplier,
                max_intron_size=args.max_intron_size,
                min_intron_size=args.min_intron_size,
                min_exon_size=args.min_exon_size,
                library_size=args.library_size,
                search_filter=args.search_filter,
                motif_search_window_size=args.motif_search_window_size,
                max_gaps_mismatches=args.max_gaps_mismatches,
                motif_radius=args.motif_radius,
                genome_bowtie1_args=args.genome_bowtie1_args,
                junction_criteria=args.junction_criteria,
                indel_criteria=args.indel_criteria,
                transcriptome_bowtie2_args=args.transcriptome_bowtie2_args,
                experimental=args.experimental,
                count_multiplier=args.count_multiplier,
                max_refs_per_strand=args.max_refs_per_strand,
                tie_margin=args.tie_margin,
                normalize_percentile=args.normalize_percentile,
                transcriptome_indexes_per_sample=\
                    args.transcriptome_indexes_per_sample,
                drop_deletions=args.drop_deletions,
                do_not_output_bam_by_chr=args.do_not_output_bam_by_chr,
                do_not_output_ave_bw_by_chr=args.do_not_output_ave_bw_by_chr,
                do_not_drop_polyA_tails=args.do_not_drop_polyA_tails,
                deliverables=args.deliverables,
                bam_basename=args.bam_basename,
                tsv_basename=args.tsv_basename,
                bed_basename=args.bed_basename, log_uri=args.log_uri,
                ami_version=args.ami_version,
                visible_to_all_users=args.visible_to_all_users,
                tags='', name=args.name,
                action_on_failure=args.action_on_failure,
                hadoop_jar=args.hadoop_jar,
                master_instance_count=args.master_instance_count,
                master_instance_type=args.master_instance_type,
                master_instance_bid_price=args.master_instance_bid_price,
                core_instance_count=args.core_instance_count,
                core_instance_type=args.core_instance_type,
                core_instance_bid_price=args.core_instance_bid_price,
                task_instance_count=args.task_instance_count,
                task_instance_type=args.task_instance_type,
                task_instance_bid_price=args.task_instance_bid_price,
                ec2_key_name=args.ec2_key_name,
                keep_alive=args.keep_alive,
                termination_protected=args.termination_protected,
                consistent_view=args.consistent_view,
                no_direct_copy=args.no_direct_copy,
                check_manifest=(not args.do_not_check_manifest),
                intermediate_lifetime=args.intermediate_lifetime,
                max_task_attempts=args.max_task_attempts,
                dbgap_key=args.dbgap_key,
                secure_stack_name=args.secure_stack_name,
                ec2_subnet_id=args.ec2_subnet_id,
                ec2_master_security_group_id=args.ec2_master_security_group_id,
                ec2_slave_security_group_id=args.ec2_slave_security_group_id
            )
    elif args.job_flow == 'align' and args.align_mode == 'elastic':
        mode = 'elastic'
        json_creator = RailRnaElasticAlignJson(
                args.manifest, args.output,
                assembly=args.assembly,
                isofrag_idx=args.isofrag_idx,
                intermediate_dir=args.intermediate,
                force=args.force, aws_exe=args.aws, profile=args.profile,
                region=args.region,
                service_role=args.service_role,
                instance_profile=args.instance_profile,
                verbose=args.verbose,
                input_dir=args.input,
                k=args.k, bowtie2_args=args.bowtie2_args,
                partition_length=args.partition_length,
                max_readlet_size=args.max_readlet_size,
                readlet_config_size=args.readlet_config_size,
                min_readlet_size=args.min_readlet_size,
                readlet_interval=args.readlet_interval,
                cap_size_multiplier=args.cap_size_multiplier,
                max_intron_size=args.max_intron_size,
                min_intron_size=args.min_intron_size,
                min_exon_size=args.min_exon_size,
                library_size=args.library_size,
                search_filter=args.search_filter,
                motif_search_window_size=args.motif_search_window_size,
                max_gaps_mismatches=args.max_gaps_mismatches,
                motif_radius=args.motif_radius,
                genome_bowtie1_args=args.genome_bowtie1_args,
                junction_criteria=args.junction_criteria,
                indel_criteria=args.indel_criteria,
                transcriptome_bowtie2_args=args.transcriptome_bowtie2_args,
                experimental=args.experimental,
                count_multiplier=args.count_multiplier,
                max_refs_per_strand=args.max_refs_per_strand,
                tie_margin=args.tie_margin,
                normalize_percentile=args.normalize_percentile,
                transcriptome_indexes_per_sample=\
                    args.transcriptome_indexes_per_sample,
                drop_deletions=args.drop_deletions,
                do_not_output_bam_by_chr=args.do_not_output_bam_by_chr,
                do_not_output_ave_bw_by_chr=args.do_not_output_ave_bw_by_chr,
                do_not_drop_polyA_tails=args.do_not_drop_polyA_tails,
                deliverables=args.deliverables,
                bam_basename=args.bam_basename, tsv_basename=args.tsv_basename,
                bed_basename=args.bed_basename, log_uri=args.log_uri,
                ami_version=args.ami_version,
                visible_to_all_users=args.visible_to_all_users,
                tags='', name=args.name,
                action_on_failure=args.action_on_failure,
                hadoop_jar=args.hadoop_jar,
                master_instance_count=args.master_instance_count,
                master_instance_type=args.master_instance_type,
                master_instance_bid_price=args.master_instance_bid_price,
                core_instance_count=args.core_instance_count,
                core_instance_type=args.core_instance_type,
                core_instance_bid_price=args.core_instance_bid_price,
                task_instance_count=args.task_instance_count,
                task_instance_type=args.task_instance_type,
                task_instance_bid_price=args.task_instance_bid_price,
                ec2_key_name=args.ec2_key_name,
                keep_alive=args.keep_alive,
                termination_protected=args.termination_protected,
                consistent_view=args.consistent_view,
                no_direct_copy=args.no_direct_copy,
                intermediate_lifetime=args.intermediate_lifetime,
                max_task_attempts=args.max_task_attempts,
                secure_stack_name=args.secure_stack_name,
                ec2_subnet_id=args.ec2_subnet_id,
                ec2_master_security_group_id=args.ec2_master_security_group_id,
                ec2_slave_security_group_id=args.ec2_slave_security_group_id
            )
    elif args.job_flow == 'prep' and args.prep_mode == 'elastic':
        mode = 'elastic'
        json_creator = RailRnaElasticPreprocessJson(
                args.manifest, args.output,
                do_not_bin_quals=args.do_not_bin_quals,
                short_read_names=args.short_read_names,
                skip_bad_records=args.skip_bad_records,
                ignore_missing_sra_samples=args.ignore_missing_sra_samples,
                intermediate_dir=args.intermediate,
                force=args.force, aws_exe=args.aws, profile=args.profile,
                region=args.region,
                service_role=args.service_role,
                instance_profile=args.instance_profile,
                verbose=args.verbose,
                log_uri=args.log_uri,
                ami_version=args.ami_version,
                visible_to_all_users=args.visible_to_all_users,
                tags='', name=args.name,
                action_on_failure=args.action_on_failure,
                hadoop_jar=args.hadoop_jar,
                master_instance_count=args.master_instance_count,
                master_instance_type=args.master_instance_type,
                master_instance_bid_price=args.master_instance_bid_price,
                core_instance_count=args.core_instance_count,
                core_instance_type=args.core_instance_type,
                core_instance_bid_price=args.core_instance_bid_price,
                task_instance_count=args.task_instance_count,
                task_instance_type=args.task_instance_type,
                task_instance_bid_price=args.task_instance_bid_price,
                ec2_key_name=args.ec2_key_name,
                keep_alive=args.keep_alive,
                termination_protected=args.termination_protected,
                consistent_view=args.consistent_view,
                no_direct_copy=args.no_direct_copy,
                check_manifest=(not args.do_not_check_manifest),
                intermediate_lifetime=args.intermediate_lifetime,
                max_task_attempts=args.max_task_attempts,
                dbgap_key=args.dbgap_key,
                secure_stack_name=args.secure_stack_name,
                ec2_subnet_id=args.ec2_subnet_id,
                ec2_master_security_group_id=args.ec2_master_security_group_id,
                ec2_slave_security_group_id=args.ec2_slave_security_group_id
            )
    # Launch
    try:
        log_file = os.path.join(
                            args.log, 
                            'flow.%s.log' % datetime.datetime.now().isoformat()
                        )
    except AttributeError:
        # No log file
        log_file = None
    try:
        region_to_use = json_creator.base.region
    except AttributeError:
        # No region specified; use US Standard
        region_to_use = 'us-east-1'
    launcher = Launcher(force=json_creator.base.force,
                                    num_processes=(
                                        json_creator.base.num_processes
                                        if mode == 'local'
                                        else None
                                    ),
                                    keep_intermediates=(
                                       args.keep_intermediates
                                       if mode in ['local', 'parallel']
                                       else False
                                    ),
                                    gzip_intermediates=(
                                       args.gzip_intermediates
                                       if mode in ['local', 'parallel']
                                       else False
                                    ),
                                    gzip_level=(
                                       args.gzip_level
                                       if mode in ['local', 'parallel']
                                       else 3
                                    ),
                                    sort_memory_cap=(
                                        args.sort_memory_cap
                                        if mode in ['local', 'parallel']
                                        else 0.2
                                    ),
                                    max_task_attempts=(
                                        args.max_task_attempts
                                        if mode in ['local', 'parallel']
                                        else 1
                                    ),
                                    log=(
                                        log_file if mode
                                        in ['local', 'parallel']
                                        else None
                                    ),
                                    region=region_to_use,
                                    common=(
                                        args.log
                                        if mode == 'parallel'
                                        else None
                                    ),
                                    ipython_profile=(
                                        args.ipython_profile
                                        if mode == 'parallel'
                                        else None
                                    ),
                                    ipcontroller_json=(
                                        args.ipcontroller_json
                                        if mode == 'parallel'
                                        else None
                                    ),
                                    scratch=(
                                        args.scratch
                                        if mode in ['local', 'parallel']
                                        else None
                                    ),
                                    direct_write=(
                                        args.direct_write
                                        if mode == 'parallel'
                                        else False
                                    ),
                                    sort=(
                                        json_creator.base.sort_exe
                                        if mode in ['local', 'parallel']
                                        else None
                                    ),
                                    json=args.json,
                                    profile=(
                                        args.profile
                                        if mode == 'elastic'
                                        else None
                                    ),
                                )
    launcher.run(mode, json.dumps(json_creator.json_serial))