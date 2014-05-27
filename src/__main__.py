#!/usr/bin/env python
"""
rail-rna
Part of Rail-RNA

Main executable for Rail-RNA. Prints introductory message and "controls"
mode arguments {local, cloud} and flow arguments {preprocess, align, all}.
Command-line interface is inspired by git's and bowtie's.
"""

import os
base_path = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
driver_path = os.path.join(base_path, 'rna', 'driver')
import site
site.addsitedir(driver_path)
site.addsitedir(base_path)
from rna_config import *
from dooplicity.tools import which
from version import version_number
import argparse
import sys
import json
import subprocess
from argparse import SUPPRESS

_usage_message = \
"""rail-rna <job flow> <mode> <[args]>

  <job flow>       {{prep, align, go}}
                     prep: preprocess reads listed in a required manifest
                       file (specified with --manifest)
                     align: align preprocessed reads (specified with --input)
                     go: perform prep and align in succession
  <mode>           {{local, elastic}}
                     local: run Rail-RNA on this computer
                     elastic: run Rail-RNA on Amazon Elastic MapReduce.
                       Requires that the user sign up for Amazon Web Services

=x= Rail-RNA v{0} by Abhi Nellore (anellore@jhu.edu; www.github.com/buci)

Rail-RNA is a scalable MapReduce pipeline that can analyze many RNA-seq
datasets at once. To view help for a given combination of <job flow> and
<mode>, specify both, then add -h/--help.""".format(version_number)
_help_set = set(['--help', '-h'])
_argv_set = set(sys.argv)

def general_usage(job_flow_and_mode, required_args=''):
    """ Special Rail-RNA usage message at subcommand level.

        job_flow_and_mode: job flow and mode separated by a space

        Return value: usage message
    """
    return \
"""rail-rna {0} {1}<[opts]>{2}""".format(
job_flow_and_mode, required_args,
"""

Add --help/-h to view help.""" if not _help_set.intersection(_argv_set)
else ''
)

class RailParser(argparse.ArgumentParser):
    """ Accommodates Rail-RNA's subcommand structure. """
    
    def error(self, message):
        if not _help_set.intersection(_argv_set):
            print >>sys.stderr, 'error: %s' % message
        self.print_usage()
        sys.exit(2)

class RailHelpFormatter(argparse.HelpFormatter):
    """ Formats help in a more condensed way.

        Overrides not-so-public argparse API, but since the Python 2.x line is
        no longer under very active development, this is probably okay.
    """

    def _get_help_string(self, action):
        help = action.help
        if '(def: ' not in action.help and not action.required \
            and not action.const:
            if action.default is not argparse.SUPPRESS:
                defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    help += ' (def: %(default)s)'
        return help

    def _format_action_invocation(self, action):
        if not action.option_strings:
            metavar, = self._metavar_formatter(action, action.dest)(1)
            return metavar
        else:
            return '%s %s' % ('/'.join(action.option_strings),
                                self._format_args(action, action.dest.upper()))

class Launcher:
    """ Facilitates replacing the current process with a Dooplicity runner. """

    def __init__(self, force=False, num_processes=1, region='us-east-1'):
        self.force = force
        self.num_processes = num_processes
        self.region = region

    def run(self, mode, payload):
        """ Replaces current process, using PyPy if it's available.

            Transferring control allows Dooplicity to handle errors from here
            on out.

            mode: 'local' or 'elastic'; launches hadoop_simulator.py or 
                emr_runner.py, respectively
            payload: string with json payload to copy to stdin
                of replacement process
        """
        read_pipe, write_pipe = os.pipe()
        if mode == 'local':
            if 'pypy 2.' in sys.version.lower():
                # Executable has the user's desired version of PyPy
                executable = sys.executable
            else:
                pypy_exe = which('pypy')
                print_warning = False
                if pypy_exe is not None:
                    try:
                        if 'pypy 2.' in \
                            subprocess.check_output(
                                    [pypy_exe, '--version'],
                                    stderr=subprocess.STDOUT
                                ).lower():
                            executable = pypy_exe
                        else:
                            executable = sys.executable
                            print_warning = True
                    except Exception:
                        executable = sys.executable
                        print_warning = True
                else:
                    print_warning = True
                if print_warning:
                    warning_message = ('WARNING: PyPy 2.x not found. '
                        'Installation is recommended to optimize performance '
                        'of Rail-RNA in "local" mode. If it is installed, '
                        'make sure the "pypy" executable is in PATH, or use '
                        'it to execute Rail-RNA via "[path to \'pypy\'] '
                        '[path to \'rail-rna\']".')
                    executable = sys.executable
                else:
                    warning_message = 'Launching Dooplicity with PyPy....'
                print >>sys.stderr, warning_message
                if not sys.stderr.isatty():
                    # So the user sees it too
                    print warning_message
            runner_args = [executable, os.path.join(base_path, 'dooplicity',
                                    'hadoop_simulator.py'),
                            '-p', str(self.num_processes),
                            '-b', os.path.join(base_path, 
                                    'rna', 'driver', 'rail-rna.txt')]
            if self.force:
                runner_args.append('-f')
        else:
            # Just run Python; who cares?
            executable = sys.executable
            runner_args = [executable, os.path.join(base_path, 'dooplicity',
                                    'emr_runner.py'),
                            '-b', os.path.join(base_path, 
                                    'rna', 'driver', 'rail-rna.txt')]
            if self.force:
                runner_args.append('-f')
            if self.region != 'us-east-1':
                runner_args.extend(['-r', self.region])
        os.write(write_pipe, payload)
        os.close(write_pipe)
        os.dup2(read_pipe, sys.stdin.fileno())
        os.execv(executable, runner_args)
        ###SCRIPT TERMINATES HERE###


def rail_help_wrapper(prog):
    """ So formatter_class's max_help_position can be changed. """
    return RailHelpFormatter(prog, max_help_position=37)

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
                                        '-m <file> -i <dir> -1 <idx> '
                                        '-2 <idx> '),
                                    formatter_class=rail_help_wrapper,
                                    add_help=False
                                )
    align_elastic_parser = align_parsers.add_parser(
                                    'elastic',
                                    usage=general_usage('align elastic',
                                        '-m <file> -i <s3_dir> -a '
                                        '<choice/tgz> \r\n       '
                                        '-c <int> -o <s3_dir> '),
                                    formatter_class=rail_help_wrapper,
                                    add_help=False
                                )
    go_parsers = go_mode_parser.add_subparsers(dest='go_mode')
    go_local_parser = go_parsers.add_parser(
                                    'local',
                                    usage=general_usage('go local',
                                        '-m <file> -1 <idx> -2 <idx> '),
                                    formatter_class=rail_help_wrapper,
                                    add_help=False
                                )
    go_elastic_parser = go_parsers.add_parser(
                                    'elastic',
                                    usage=general_usage('go elastic',
                                        '-m <file> -i <s3_dir> -a '
                                        '<choice/s3_dir> \r\n       '
                                        '-o <s3_dir> '),
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
                        go_local_general, prep_elastic_general,
                        align_elastic_general, go_elastic_general]:
        subparser.add_argument(
                    '-h', '--help',
                    action='help', default=SUPPRESS,
                    help='show this help message and exit'
                )
    RailRnaErrors.add_args(general_parser=go_elastic_general,
                            exec_parser=go_elastic_exec,
                            required_parser=go_elastic_required)
    RailRnaErrors.add_args(general_parser=go_local_general,
                            exec_parser=go_local_exec,
                            required_parser=go_local_required),
    RailRnaErrors.add_args(general_parser=align_elastic_general,
                            exec_parser=align_elastic_exec,
                            required_parser=align_elastic_required)
    RailRnaErrors.add_args(general_parser=align_local_general,
                            exec_parser=align_local_exec,
                            required_parser=align_local_required)
    RailRnaErrors.add_args(general_parser=prep_elastic_general,
                            exec_parser=prep_elastic_exec,
                            required_parser=prep_elastic_required)
    RailRnaErrors.add_args(general_parser=prep_local_general,
                            exec_parser=prep_local_exec,
                            required_parser=prep_local_required)
    RailRnaLocal.add_args(required_parser=go_local_required,
                            general_parser=go_local_general,
                            output_parser=go_local_output, 
                            prep=False, align=False)
    RailRnaLocal.add_args(required_parser=align_local_required,
                            general_parser=align_local_general,
                            output_parser=align_local_output,
                            prep=False, align=True)
    RailRnaLocal.add_args(required_parser=prep_local_required,
                            general_parser=prep_local_general,
                            output_parser=prep_local_output,
                            prep=True, align=False)
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
                               output_parser=prep_elastic_output)
    RailRnaPreprocess.add_args(general_parser=prep_local_general,
                               output_parser=prep_local_output)
    RailRnaPreprocess.add_args(general_parser=go_elastic_general,
                               output_parser=go_elastic_general)
    RailRnaPreprocess.add_args(general_parser=go_local_general,
                               output_parser=go_local_general)
    RailRnaAlign.add_args(required_parser=align_local_required,
                          exec_parser=align_local_exec,
                          output_parser=align_local_output,
                          algo_parser=align_local_algo, elastic=False)
    RailRnaAlign.add_args(required_parser=align_elastic_required,
                          exec_parser=align_elastic_exec,
                          output_parser=align_elastic_output,
                          algo_parser=align_elastic_algo, elastic=True)
    RailRnaAlign.add_args(required_parser=go_local_required,
                          exec_parser=go_local_exec,
                          output_parser=go_local_output,
                          algo_parser=go_local_algo, elastic=False)
    RailRnaAlign.add_args(required_parser=go_elastic_required,
                          exec_parser=go_elastic_exec,
                          output_parser=go_elastic_output,
                          algo_parser=go_elastic_algo, elastic=True)
    args = parser.parse_args()
    if args.job_flow == 'go' and args.go_mode == 'local':
        mode = 'local'
        json_creator = RailRnaLocalAllJson(
                args.manifest, args.output,
                intermediate_dir=args.intermediate,
                force=args.force, aws_exe=args.aws, profile=args.profile,
                verbose=args.verbose,
                nucleotides_per_input=args.nucleotides_per_input,
                gzip_input=(not args.do_not_gzip_input),
                bowtie1_idx=args.bowtie1_idx, bowtie2_idx=args.bowtie2_idx,
                bowtie1_exe=args.bowtie1, bowtie2_exe=args.bowtie2,
                bowtie1_build_exe=args.bowtie1_build,
                bowtie2_build_exe=args.bowtie2_build,
                bowtie2_args=args.bowtie2_args,
                samtools_exe=args.samtools,
                bedtobigbed_exe=args.bedtobigbed,
                genome_partition_length=args.genome_partition_length,
                max_readlet_size=args.max_readlet_size,
                min_readlet_size=args.min_readlet_size,
                readlet_interval=args.readlet_interval,
                cap_size_multiplier=args.cap_size_multiplier,
                max_intron_size=args.max_intron_size,
                min_intron_size=args.min_intron_size,
                min_exon_size=args.min_exon_size,
                motif_search_window_size=args.motif_search_window_size,
                motif_radius=args.motif_radius,
                normalize_percentile=args.normalize_percentile,
                do_not_output_bam_by_chr=args.do_not_output_bam_by_chr,
                output_sam=args.output_sam, bam_basename=args.bam_basename,
                bed_basename=args.bed_basename,
                num_processes=args.num_processes,
                keep_intermediates=args.keep_intermediates,
                check_manifest=(not args.do_not_check_manifest)
            )
    elif args.job_flow == 'align' and args.align_mode == 'local':
        mode = 'local'
        json_creator = RailRnaLocalAlignJson(
                args.manifest, args.output, args.input,
                intermediate_dir=args.intermediate,
                force=args.force, aws_exe=args.aws, profile=args.profile,
                verbose=args.verbose,
                bowtie1_idx=args.bowtie1_idx, bowtie2_idx=args.bowtie2_idx,
                bowtie1_exe=args.bowtie1, bowtie2_exe=args.bowtie2,
                bowtie1_build_exe=args.bowtie1_build,
                bowtie2_build_exe=args.bowtie2_build,
                bowtie2_args=args.bowtie2_args,
                samtools_exe=args.samtools,
                bedtobigbed_exe=args.bedtobigbed,
                genome_partition_length=args.genome_partition_length,
                max_readlet_size=args.max_readlet_size,
                min_readlet_size=args.min_readlet_size,
                readlet_interval=args.readlet_interval,
                cap_size_multiplier=args.cap_size_multiplier,
                max_intron_size=args.max_intron_size,
                min_intron_size=args.min_intron_size,
                min_exon_size=args.min_exon_size,
                motif_search_window_size=args.motif_search_window_size,
                motif_radius=args.motif_radius,
                normalize_percentile=args.normalize_percentile,
                do_not_output_bam_by_chr=args.do_not_output_bam_by_chr,
                output_sam=args.output_sam, bam_basename=args.bam_basename,
                bed_basename=args.bed_basename,
                num_processes=args.num_processes,
                keep_intermediates=args.keep_intermediates
            )
    elif args.job_flow == 'prep' and args.prep_mode == 'local':
        mode = 'local'
        json_creator = RailRnaLocalPreprocessJson(
                args.manifest, args.output,
                intermediate_dir=args.intermediate,
                force=args.force, aws_exe=args.aws, profile=args.profile,
                verbose=args.verbose,
                nucleotides_per_input=args.nucleotides_per_input,
                gzip_input=(not args.do_not_gzip_input),
                num_processes=args.num_processes,
                keep_intermediates=args.keep_intermediates
            )
    elif args.job_flow == 'go' and args.go_mode == 'elastic':
        mode = 'elastic'
        json_creator = RailRnaElasticAllJson(
                args.manifest, args.output,
                intermediate_dir=args.intermediate,
                force=args.force, aws_exe=args.aws, profile=args.profile,
                region=args.region, verbose=args.verbose,
                nucleotides_per_input=args.nucleotides_per_input,
                gzip_input=(not args.do_not_gzip_input),
                bowtie2_args=args.bowtie2_args,
                genome_partition_length=args.genome_partition_length,
                max_readlet_size=args.max_readlet_size,
                min_readlet_size=args.min_readlet_size,
                readlet_interval=args.readlet_interval,
                cap_size_multiplier=args.cap_size_multiplier,
                max_intron_size=args.max_intron_size,
                min_intron_size=args.min_intron_size,
                min_exon_size=args.min_exon_size,
                motif_search_window_size=args.motif_search_window_size,
                motif_radius=args.motif_radius,
                normalize_percentile=args.normalize_percentile,
                do_not_output_bam_by_chr=args.do_not_output_bam_by_chr,
                output_sam=args.output_sam, bam_basename=args.bam_basename,
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
                check_manifest=(not args.do_not_check_manifest)
            )
    elif args.job_flow == 'align' and args.align_mode == 'elastic':
        mode = 'elastic'
        json_creator = RailRnaElasticAlignJson(
                args.manifest, args.output,
                intermediate_dir=args.intermediate,
                force=args.force, aws_exe=args.aws, profile=args.profile,
                region=args.region, verbose=args.verbose,
                input_dir=args.input,
                bowtie2_args=args.bowtie2_args,
                genome_partition_length=args.genome_partition_length,
                max_readlet_size=args.max_readlet_size,
                min_readlet_size=args.min_readlet_size,
                readlet_interval=args.readlet_interval,
                cap_size_multiplier=args.cap_size_multiplier,
                max_intron_size=args.max_intron_size,
                min_intron_size=args.min_intron_size,
                min_exon_size=args.min_exon_size,
                motif_search_window_size=args.motif_search_window_size,
                motif_radius=args.motif_radius,
                normalize_percentile=args.normalize_percentile,
                do_not_output_bam_by_chr=args.do_not_output_bam_by_chr,
                output_sam=args.output_sam, bam_basename=args.bam_basename,
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
                termination_protected=args.termination_protected
            )
    elif args.job_flow == 'prep' and args.prep_mode == 'elastic':
        mode = 'elastic'
        json_creator = RailRnaElasticPreprocessJson(
                args.manifest, args.output,
                intermediate_dir=args.intermediate,
                force=args.force, aws_exe=args.aws, profile=args.profile,
                region=args.region, verbose=args.verbose,
                nucleotides_per_input=args.nucleotides_per_input,
                gzip_input=(not args.do_not_gzip_input),
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
                check_manifest=(not args.do_not_check_manifest)
            )
    # Launch
    try:
        launcher = Launcher(force=json_creator.base.force,
                                        num_processes=(
                                            json_creator.base.num_processes
                                            if mode == 'local'
                                            else None
                                        ),
                                        region=args.region)
    except AttributeError:
        launcher = Launcher(force=json_creator.base.force,
                                        num_processes=(
                                            json_creator.base.num_processes
                                            if mode == 'local'
                                            else None
                                        ))
    launcher.run(mode, json.dumps(json_creator.json_serial))