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
from dooplicity.tools import path_join, is_exe, which, register_cleanup, \
    apply_async_with_errors, engine_string_from_list, cd
from version import version_number
import sys
import argparse
import subprocess
import time
from traceback import format_exc
from collections import defaultdict
import random
import string
import socket
import exe_paths
import unicodedata
from distutils.util import strtobool

_help_set = set(['--help', '-h'])
_argv_set = set(sys.argv)
_whitespace_and_comma = string.whitespace + ','

class RailParser(argparse.ArgumentParser):
    """ Accommodates Rail-RNA's subcommand structure. """
    
    def error(self, message):
        if not _help_set.intersection(_argv_set):
            print >>sys.stderr, 'error: %s' % message
            self.print_usage()
            sys.exit(2)
        self.print_usage()
        sys.exit(0)

class RailHelpFormatter(argparse.HelpFormatter):
    """ Formats help in a more condensed way.

        Overrides not-so-public argparse API, but since the Python 2.x line is
        no longer under very active development, this is probably okay.
    """

    def _get_help_string(self, action):
        action_help = action.help
        if '(def: ' not in action.help and not action.required \
            and not action.const:
            if action.default is not argparse.SUPPRESS:
                defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE,
                                    argparse.ONE_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    action_help += ' (def: %(default)s)'
        # Take out the "def: "
        return action_help

    def _format_action_invocation(self, action):
        if not action.option_strings:
            metavar, = self._metavar_formatter(action, action.dest)(1)
            return metavar
        else:
            return '%s %s' % ('/'.join(action.option_strings),
                                self._format_args(action, action.dest.upper()))

    def _format_args(self, action, default_metavar):
        get_metavar = self._metavar_formatter(action, default_metavar)
        if action.nargs in [None, argparse.ZERO_OR_MORE, argparse.ONE_OR_MORE]:
            result = '%s' % get_metavar(1)
        elif action.nargs == argparse.OPTIONAL:
            result = '[%s]' % get_metavar(1)
        elif action.nargs == argparse.REMAINDER:
            result = '...'
        elif action.nargs == argparse.PARSER:
            result = '%s ...' % get_metavar(1)
        else:
            formats = ['%s' for _ in range(action.nargs)]
            result = ' '.join(formats) % get_metavar(action.nargs)
        return result

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

def rail_help_wrapper(prog):
    """ So formatter_class's max_help_position can be changed. """
    return RailHelpFormatter(prog, max_help_position=40)

'''These are placed here for convenience; their locations may change
on EMR depending on bootstraps.'''
_hadoop_streaming_jar = '/home/hadoop/contrib/streaming/hadoop-streaming.jar'
_jar_target = '/mnt/lib'
_multiple_files_jar = _jar_target + '/multiple-files.jar'
_relevant_elephant_jar = _jar_target + '/relevant-elephant.jar'
_mod_partitioner_jar = _jar_target + '/mod-partitioner.jar'
_hadoop_lzo_jar = ('/home/hadoop/.versions/2.4.0/share/hadoop'
                   '/common/lib/hadoop-lzo.jar')
_s3distcp_jar = '/home/hadoop/lib/emr-s3distcp-1.0.jar'
_hdfs_temp_dir = 'hdfs:///railtemp'
_base_combine_split_size = 268435456 # 250 MB
_elastic_bowtie1_idx = '/mnt/space/index/genome'
_elastic_bowtie2_idx = '/mnt/space/index/genome'
_elastic_bedgraphtobigwig_exe = 'bedGraphToBigWig'
_elastic_samtools_exe = 'samtools'
_elastic_bowtie1_exe = 'bowtie'
_elastic_bowtie2_exe = 'bowtie2'
_elastic_bowtie1_build_exe = 'bowtie-build'
_elastic_bowtie2_build_exe = 'bowtie2-build'
_elastic_fastq_dump_exe = 'fastq-dump'
_elastic_vdb_config_exe = '/usr/local/bin/vdb-config'
_elastic_vdb_workspace = '/mnt/space/sra_workspace'
_elastic_step_dir = '/usr/local/raildotbio/rail-rna/rna/steps'

_assemblies = {
        'hg19' : 's3://rail-emr-requester-pays/index/hg19.tar.gz',
        'hg38' : 's3://rail-emr-requester-pays/index/hg38.tar.gz',
        'mm9' : 's3://rail-emr-requester-pays/index/mm9.tar.gz',
        'mm10' : 's3://rail-emr-requester-pays/index/mm10.tar.gz',
        'dm3' : 's3://rail-emr-requester-pays/index/dm3.tar.gz',
        'dm6' : 's3://rail-emr-requester-pays/index/dm6.tar.gz'
    }

# Set basename of the transcript fragment index; can't settle on this
_transcript_fragment_idx_basename = 'isofrags'

# Decide Python executable
if 'pypy 2.' in sys.version.lower():
    # Executable has the user's desired version of PyPy
    _warning_message = 'Launching Dooplicity runner with PyPy...'
    _executable = sys.executable
else:
    if exe_paths.pypy is None:
        _pypy_exe = which('pypy')
    else:
        candidate_pypy_exe = which(exe_paths.pypy)
        if candidate_pypy_exe is not None:
            _pypy_exe = exe_paths.pypy
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

def print_to_screen(message, newline=True, carriage_return=False):
    """ Prints message to stdout as well as stderr if stderr is redirected.

        message: message to print
        newline: True iff newline should be printed
        carriage_return: True iff carriage return should be printed; also
            clears line with ANSI escape code

        No return value.
    """
    full_message = ('\x1b[K' + message + ('\r' if carriage_return else '')
                        + ('\n' if newline else ''))
    try:
        sys.stderr.write(full_message)
        if sys.stderr.isatty():
            sys.stderr.flush()
        else:
            try:
                # So the user sees it too
                sys.stdout.write(full_message)
                sys.stdout.flush()
            except UnicodeEncodeError:
                sys.stdout.write(
                                unicodedata.normalize(
                                        'NFKD', full_message
                                    ).encode('ascii', 'ignore')
                            )
                sys.stdout.flush()
    except UnicodeEncodeError:
        sys.stderr.write(
                        unicodedata.normalize(
                                'NFKD', full_message
                            ).encode('ascii', 'ignore')
                    )
        sys.stderr.flush()

def ready_engines(rc, base, prep=False):
    """ Prepares engines for checks and copies Rail/manifest/index to nodes. 

        rc: IPython Client object
        base: instance of RailRnaErrors
        prep: True iff it's a preprocess job flow

        No return value.
    """
    try:
        import IPython
    except ImportError:
        # Should have been taken care of by a different fxn, but just in case
        raise RuntimeError(
               'IPython is required to run Rail-RNA in '
               '"parallel" mode. Visit ipython.org to '
               'download it, or simply download the Anaconda '
               'distribution of Python at '
               'https://store.continuum.io/cshop/anaconda/; it\'s '
               'easy to install and comes with IPython and '
               'several other useful packages.'
            )
    all_engines = rc.ids
    '''Clear remote namespaces.'''
    rc[:].clear()
    '''Test that intermediate directory is accessible from everywhere; create
    dir in process.'''
    try:
        os.makedirs(base.intermediate_dir)
    except OSError:
        # Hopefully exists
        pass
    # Create dud file in intermediate directory
    dud_filename = os.path.join(base.intermediate_dir, 
                            ''.join(random.choice(string.ascii_uppercase
                                    + string.digits) for _ in xrange(40))
                        )
    with open(dud_filename, 'w') as dud_stream:
        print >>dud_stream, 'DUD'
    '''Now test for existence of dud file across engines; a dud file with a
    very random name is created to ensure that the directory being searched for
    is really the one specified by the user rather than some directory that 
    only an engine could see, which with high probability is absent this
    file.'''
    try:
        dud_results = apply_async_with_errors(rc, all_engines, os.path.exists,
                                            dud_filename, dict_format=True,
                                            message=('Error(s) encountered '
                                                     'testing that '
                                                     'the log directory is '
                                                     'accessible from '
                                                     'all engines. Restart '
                                                     'IPython engines '
                                                     'and try again.')
                                        )
    finally:
        # No matter what, kill the dud
        os.remove(dud_filename)
    bad_engines = [engine for engine in dud_results if not dud_results[engine]]
    if bad_engines:
        raise RuntimeError(('Engines %s cannot access the log directory %s. '
                            'Ensure that the log directory is in a location '
                            'accessible from all engines.') % (
                                    engine_string_from_list(bad_engines),
                                    base.intermediate_dir
                                ))
    current_hostname = socket.gethostname()
    engine_to_hostnames = apply_async_with_errors(
                                rc, all_engines, socket.gethostname,
                                dict_format=True
                            )
    hostname_to_engines = defaultdict(set)
    for engine in engine_to_hostnames:
        hostname_to_engines[engine_to_hostnames[engine]].add(engine)
    '''Select engines to do "heavy lifting"; that is, they remove files copied
    to hosts on SIGINT/SIGTERM. Do it randomly (NO SEED) so if IWF occurs,
    second try will be different. IWF = intermittent weird failure, terminology 
    borrowed from a PC repair guide from the nineties that one of us (AN) wants
    to perpetuate.'''
    pids = apply_async_with_errors(rc, all_engines, os.getpid)
    # Set random seed so temp directory is reused if restarting Rail
    random.seed(str(sorted(pids)))
    engines_for_copying = [random.choice(list(engines)) 
                            for engines in hostname_to_engines.values()
                            if len(engines) > 0]
    '''Herd won't work with local engines, work around this by separating
    engines into two groups: local and remote.'''
    remote_hostnames_for_copying = list(
            set(hostname_to_engines.keys()).difference(set([current_hostname]))
        )
    local_engines_for_copying = [engine for engine in engines_for_copying
                                 if engine
                                 in hostname_to_engines[current_hostname]]
    '''Create temporary directories on selected nodes; NOT WINDOWS-COMPATIBLE;
    must be changed if porting Rail to Windows.'''
    if base.scratch is None:
        scratch_dir = tempfile.gettempdir()
    else:
        scratch_dir = base.scratch
    temp_dir = os.path.join(scratch_dir, 'railrna-%s' %
                            ''.join(random.choice(string.ascii_uppercase
                                    + string.digits) for _ in xrange(12)))
    if not prep and not base.do_not_copy_index_to_nodes:
        dir_to_create = os.path.join(temp_dir, 'genome')
    else:
        dir_to_create = temp_dir
    temp_dirs = apply_async_with_errors(rc, all_engines,
        subprocess.check_output,
        'echo "%s"' % temp_dir,
        shell=True,
        executable='/bin/bash',
        message=('Error obtaining full paths of temporary directories '
                 'on cluster nodes. Restart IPython engines '
                 'and try again.'
            ),
        dict_format=True)
    for engine in temp_dirs:
        temp_dirs[engine] = temp_dirs[engine].strip()
    engines_with_unique_scratch, engines_to_symlink = [], []
    engine_to_copy_engine = {}
    for engine_for_copying in engines_for_copying:
        for engine in hostname_to_engines[
                    engine_to_hostnames[engine_for_copying]
                ]:
            engine_to_copy_engine[engine] = engine_for_copying
            if (engine != engine_for_copying
                and temp_dirs[engine] != temp_dirs[engine_for_copying]):
                engines_with_unique_scratch.append(engine)
                engines_to_symlink.append(engine)
            elif engine == engine_for_copying:
                engines_with_unique_scratch.append(engine)
    '''Must use mkdir -p rather than os.makedirs so node-specific BASH
    variable, if it was specified, works.'''
    apply_async_with_errors(rc, engines_for_copying,
        subprocess.check_output,
        'mkdir -p %s' % dir_to_create, shell=True, executable='/bin/bash',
        message=('Error(s) encountered creating temporary '
                 'directories for storing Rail on slave nodes. '
                 'Restart IPython engines and try again.'))
    if engines_to_symlink:
        # Create symlinks to resources in case of slot-local scratch dirs
        print_to_screen(
                'Adding symlinks...',
                newline=False, carriage_return=True
            )
        source_paths, destination_paths = {}, {}
        for engine_to_symlink in engines_to_symlink:
            source_paths[engine_to_symlink] = temp_dirs[
                        engine_to_copy_engine[engine_to_symlink]
                    ]
            destination_paths[engine_to_symlink] = temp_dirs[engine_to_symlink]
        apply_async_with_errors(rc, engines_to_symlink,
            os.remove, destination_paths,
            message=('Error(s) encountered removing symlinks '
                     'in slot-local scratch directories.'),
            errors_to_ignore=['OSError'])
        apply_async_with_errors(rc, engines_to_symlink,
            os.symlink, source_paths, destination_paths,
            message=('Error(s) encountered symlinking '
                     'among slot-local scratch directories.'))
        print_to_screen(
                'Added symlinks.',
                newline=True, carriage_return=True
            )
    '''Only foolproof way to die is by process polling. See
    http://stackoverflow.com/questions/284325/
    how-to-make-child-process-die-after-parent-exits for more information.'''
    apply_async_with_errors(rc, engines_with_unique_scratch,
        subprocess.check_output,
        ('echo "trap \\"{{ rm -rf {temp_dir}; exit 0; }}\\" '
         'SIGHUP SIGINT SIGTERM EXIT; '
         'while [[ \$(ps -p \$\$ -o ppid=) -gt 1 ]]; do sleep 1; done & wait" '
         '>{temp_dir}/delscript.sh').format(temp_dir=temp_dir),
        shell=True,
        executable='/bin/bash',
        message=(
                'Error creating script for scheduling temporary directories '
                'on cluster nodes for deletion. Restart IPython engines '
                'and try again.'
            ))
    apply_async_with_errors(rc, engines_with_unique_scratch, subprocess.Popen,
            '/usr/bin/env bash %s/delscript.sh' % temp_dir, shell=True,
            executable='/bin/bash',
            message=(
                'Error scheduling temporary directories on slave nodes '
                'for deletion. Restart IPython engines and try again.'
            ))
    # Compress Rail-RNA and distribute it to nodes
    compressed_rail_file = 'rail.tar.gz'
    compressed_rail_path = os.path.join(os.path.abspath(base.intermediate_dir),
                                            compressed_rail_file)
    compressed_rail_destination = os.path.join(temp_dir, compressed_rail_file)
    import tarfile
    with tarfile.open(compressed_rail_path, 'w:gz') as tar_stream:
        tar_stream.add(base_path, arcname='rail')
    try:
        import herd.herd as herd
        if '$' in compressed_rail_destination: raise ImportError
    except ImportError:
        '''Torrent distribution channel for compressed archive not available,
        or we need to use slot-local BASH variables.'''
        print_to_screen('Copying Rail-RNA to cluster nodes...',
                            newline=False, carriage_return=True)
        apply_async_with_errors(rc, engines_for_copying,
            subprocess.check_output, 'cp %s %s' % (
                    compressed_rail_path, compressed_rail_destination
                ), shell=True, executable='/bin/bash',
            message=('Error(s) encountered copying Rail to '
                     'slave nodes. Refer to the errors above -- and '
                     'especially make sure "%s" is not out of space on any '
                     'node supporting an IPython engine '
                     '-- before trying again.') % temp_dir,
        )
        print_to_screen('Copied Rail-RNA to cluster nodes.',
                            newline=True, carriage_return=False)
    else:
        if local_engines_for_copying:
            print_to_screen('Copying Rail-RNA to local filesystem...',
                            newline=False, carriage_return=True)
            apply_async_with_errors(rc, local_engines_for_copying,
                shutil.copyfile, compressed_rail_path,
                compressed_rail_destination,
                message=('Error(s) encountered copying Rail to '
                         'local filesystem. Refer to the errors above -- and '
                         'especially make sure "%s" is not out of space on '
                         'any node supporting an IPython engine '
                         '-- before trying again.') % temp_dir,
            )
            print_to_screen('Copied Rail-RNA to local filesystem.',
                                newline=True, carriage_return=False)
        if remote_hostnames_for_copying:
            print_to_screen('Copying Rail-RNA to remote nodes with Herd...')
            herd.run_with_opts(
                    compressed_rail_path,
                    compressed_rail_destination,
                    hostlist=','.join(remote_hostnames_for_copying)
                )
            print_to_screen('Copied Rail-RNA to remote nodes with Herd.',
                                newline=True, carriage_return=False)
    # Extract Rail
    print_to_screen('Extracting Rail-RNA on cluster nodes...',
                            newline=False, carriage_return=True)
    apply_async_with_errors(rc, engines_for_copying, subprocess.check_output,
            'tar xzf {} -C {}'.format(compressed_rail_destination, temp_dir),
            shell=True, executable='/bin/bash')
    print_to_screen('Extracted Rail-RNA on cluster nodes.',
                            newline=True, carriage_return=False)
    '''Add Rail to path on every engine. Must accommodate potentially different
    paths on different engines.'''
    temp_base_paths, temp_rna_paths, \
        temp_utils_paths, temp_driver_paths = {}, {}, {}, {}
    for engine in temp_dirs:
        temp_base_paths[engine] = os.path.join(temp_dirs[engine], 'rail')
        temp_utils_paths[engine] = os.path.join(temp_base_paths[engine],
                                                'rna', 'utils')
        temp_driver_paths[engine] = os.path.join(temp_base_paths[engine],
                                                 'rna', 'driver')
    apply_async_with_errors(rc, all_engines, site.addsitedir, temp_base_paths)
    apply_async_with_errors(rc, all_engines, site.addsitedir, temp_utils_paths)
    apply_async_with_errors(rc, all_engines, site.addsitedir,
                                temp_driver_paths)
    # Change to current path on every engine
    current_path = os.path.abspath('./')
    apply_async_with_errors(rc, all_engines, os.chdir, current_path,
                                errors_to_ignore=['OSError'])
    # Copy manifest to nodes
    manifest_destination = os.path.join(temp_dir, 'MANIFEST')
    try:
        import herd.herd as herd
        if '$' in manifest_destination: raise ImportError
    except ImportError:
        print_to_screen('Copying file manifest to cluster nodes...',
                            newline=False, carriage_return=True)
        apply_async_with_errors(rc, engines_for_copying, 
            subprocess.check_output, 'cp %s %s' % (
                    base.manifest, manifest_destination,
                ), shell=True, executable='/bin/bash',
            message=('Error(s) encountered copying manifest to '
                     'slave nodes. Refer to the errors above -- and '
                     'especially make sure "%s" is not out of space on any '
                     'node supporting an IPython engine '
                     '-- before trying again.') % temp_dir,
        )
        print_to_screen('Copied file manifest to cluster nodes.',
                            newline=True, carriage_return=False)
    else:
        if local_engines_for_copying:
            print_to_screen('Copying manifest to local filesystem...',
                            newline=False, carriage_return=True)
            apply_async_with_errors(rc, local_engines_for_copying,
                shutil.copyfile, base.manifest, manifest_destination,
                message=('Error(s) encountered copying manifest to '
                         'slave nodes. Refer to the errors above -- and '
                         'especially make sure "%s" is not out of space on '
                         'any node supporting an IPython engine '
                         '-- before trying again.') % temp_dir,
            )
            print_to_screen('Copied manifest to local filesystem.',
                                newline=True, carriage_return=False)
        if remote_hostnames_for_copying:
            print_to_screen('Copying manifest to remote nodes with Herd...')
            herd.run_with_opts(
                    base.manifest,
                    manifest_destination,
                    hostlist=','.join(remote_hostnames_for_copying)
                )
            print_to_screen('Copied manifest to remote nodes with Herd.',
                                newline=True, carriage_return=False)
    base.old_manifest = base.manifest
    base.manifest = manifest_destination
    if not prep and not base.do_not_copy_index_to_nodes:
        index_files = ([base.bowtie2_idx + extension
                        for extension in ['.1.bt2', '.2.bt2',
                                          '.3.bt2', '.4.bt2', 
                                          '.rev.1.bt2', '.rev.2.bt2']]
                        + [base.bowtie1_idx + extension
                            for extension in [
                                    '.1.ebwt', '.2.ebwt', '.3.ebwt',
                                    '.4.ebwt', '.rev.1.ebwt', '.rev.2.ebwt'
                                ]])
        try:
            import herd.herd as herd
            if '$' in temp_dir: raise ImportError
        except ImportError:
            print_to_screen('Warning: Herd is not installed or BASH variables '
                            'are in --scratch, so copying '
                            'Bowtie indexes to cluster nodes may be slow. '
                            'Install Herd to enable torrent distribution of '
                            'indexes across nodes, or invoke '
                            '--do-not-copy-index-to-nodes to avoid copying '
                            'indexes, which may then slow down alignment.',
                            newline=True, carriage_return=False)
            files_copied = 0
            print_to_screen(
                    'Copying Bowtie index files to cluster nodes '
                    '(%d/%d files copied)...'
                    % (files_copied, len(index_files)),
                    newline=False, carriage_return=True
                )
            for index_file in index_files:
                apply_async_with_errors(rc, engines_for_copying,
                    subprocess.check_output, 'cp %s %s' % (
                        os.path.abspath(index_file),
                        os.path.join(temp_dir, 'genome',
                                        os.path.basename(index_file))
                    ), shell=True, executable='/bin/bash',
                    message=('Error(s) encountered copying Bowtie indexes to '
                             'cluster nodes. Refer to the errors above -- and '
                             'especially make sure "%s" is not out of '
                             'space on any node supporting an IPython engine '
                             '-- before trying again.') % temp_dir
                )
                files_copied += 1
                print_to_screen(
                    'Copying Bowtie index files to cluster nodes '
                    '(%d/%d files copied)...'
                    % (files_copied, len(index_files)),
                    newline=False, carriage_return=True
                )
            print_to_screen('Copied Bowtie indexes to cluster nodes.',
                                newline=True, carriage_return=False)
        else:
            if local_engines_for_copying:
                files_copied = 0
                print_to_screen('Copying Bowtie indexes to local '
                                'filesystem (%d/%d files copied)...'
                                % (files_copied, len(index_files)),
                                newline=False, carriage_return=True)
                for index_file in index_files:
                    apply_async_with_errors(rc, engines_for_copying,
                        shutil.copyfile, os.path.abspath(index_file),
                        os.path.join(temp_dir, 'genome',
                                        os.path.basename(index_file)),
                        message=('Error(s) encountered copying Bowtie '
                                 'indexes to local filesystem. Refer to the '
                                 'errors above -- and especially make sure '
                                 '"%s" is not out of space '
                                 'on any node supporting an IPython engine '
                                 '-- before trying again.') % temp_dir
                    )
                    files_copied += 1
                    print_to_screen(
                        'Copying Bowtie indexes to local filesystem '
                        '(%d/%d files copied)...'
                        % (files_copied, len(index_files)),
                        newline=False, carriage_return=True
                    )
                print_to_screen('Copied Bowtie indexes to local '
                                'filesystem.',
                                newline=True, carriage_return=False)
            if remote_hostnames_for_copying:
                print_to_screen('Copying Bowtie indexes to cluster nodes '
                                'with Herd...',
                                newline=False, carriage_return=True)
                try:
                    for index_file in index_files:
                        herd.run_with_options(
                                os.path.abspath(index_file),
                                os.path.join(temp_dir, 'genome',
                                    os.path.basename(index_file)),
                                hostlist=','.join(hostname_to_engines.keys())
                            )
                except:
                    print_to_screen('Herd copy failed; copying without Herd.',
                                        newline=True, carriage_return=False)
                    print_to_screen('Copying Bowtie indexes to local '
                                'filesystem (%d/%d files copied)...'
                                % (files_copied, len(index_files)),
                                newline=False, carriage_return=True)
                    for index_file in index_files:
                        apply_async_with_errors(rc, engines_for_copying,
                            subprocess.check_output, 'cp %s %s' % (
                                os.path.abspath(index_file),
                                os.path.join(temp_dir, 'genome',
                                                os.path.basename(index_file))
                            ), shell=True, executable='/bin/bash',
                            message=('Error(s) encountered copying Bowtie '
                                     'indexes to local filesystem. Refer to '
                                     'the errors above -- and especially make '
                                     'sure "%s" is not out of space '
                                     'on any node supporting an IPython '
                                     'engine -- before trying again.')
                                        % temp_dir
                        )
                        files_copied += 1
                        print_to_screen(
                            'Copying Bowtie indexes to local filesystem '
                            '(%d/%d files copied)...'
                            % (files_copied, len(index_files)),
                            newline=False, carriage_return=True
                        )
                    print_to_screen('Copied Bowtie indexes to local '
                                    'filesystem.',
                                    newline=True, carriage_return=False)
                else:
                    print_to_screen('Copied Bowtie indexes to cluster nodes '
                                    'with Herd.',
                                    newline=True, carriage_return=False)
        base.bowtie1_idx = os.path.join(temp_dir, 'genome',
                                        os.path.basename(base.bowtie1_idx))
        base.bowtie2_idx = os.path.join(temp_dir, 'genome',
                                        os.path.basename(base.bowtie2_idx))
    return os.path.join(temp_dir, 'rail')

def step(name, inputs, output,
    mapper='org.apache.hadoop.mapred.lib.IdentityMapper',
    reducer='org.apache.hadoop.mapred.lib.IdentityReducer', 
    action_on_failure='TERMINATE_JOB_FLOW', jar=_hadoop_streaming_jar,
    tasks=0, partition_options=None, sort_options=None, archives=None,
    files=None, multiple_outputs=False, mod_partitioner=False,
    inputformat=None, extra_args=[]):
    """ Outputs JSON for a given step.

        name: name of step
        inputs: list of input directories/files
        output: output directory
        mapper: mapper command
        reducer: reducer command
        jar: path to Hadoop Streaming jar; ignored in local mode
        tasks: reduce task count
        partition_options: sort-like partition options or None if mapper
        sort_options: UNIX sort options or None if mapper
        archives: -archives option
        files: -files option
        multiple_outputs: True iff there are multiple outputs; else False
        mod_partitioner: True iff the mod partitioner should be used for
            the step; this partitioner assumes the key is a tuple of integers
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
    # Get key fields as max of all numbers
    if partition_options is not None:
        if sort_options is None:
            sort_options = ''
        # number of key fields is max of all numbers passed to both partitioner
        # and sorter
        key_fields = max([int(el) for arg in (partition_options.split('-k')
                                              + sort_options.split('-k'))
                            if arg.strip() and len(arg.split(',')) <= 2
                            for el in arg.strip().strip('nr').split(',')])
        to_return['HadoopJarStep']['Args'].extend([
            '-D', 'stream.num.map.output.key.fields=%d' % key_fields,
            '-D', 'mapreduce.partition.keypartitioner.options=%s'
                        % partition_options
        ])
        if sort_options:
            to_return['HadoopJarStep']['Args'].extend([
                '-D', 'mapreduce.job.output.key.comparator.class='
                      'org.apache.hadoop.mapred.lib.KeyFieldBasedComparator',
                '-D', 'mapreduce.partition.keycomparator.options=%s' % (
                                                    sort_options
                                                )
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
    if mod_partitioner:
        to_return['HadoopJarStep']['Args'][-1] \
            +=  (',%s' % _mod_partitioner_jar)
    if archives is not None:
        to_return['HadoopJarStep']['Args'].extend([
                '-archives', archives
            ])
    if files is not None:
        to_return['HadoopJarStep']['Args'].extend([
                '-files', files
            ])
    if mod_partitioner:
        to_return['HadoopJarStep']['Args'].extend([
                '-partitioner',
                'edu.jhu.cs.ModPartitioner',
            ])
    else:
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
def steps(protosteps, action_on_failure, jar, step_dir, reducer_count,
            intermediate_dir, extra_args=[], unix=False, no_direct_copy=False):
    """ Turns list with "protosteps" into well-formed StepConfig list.

        A protostep looks like this:
            {
                'name' : [name of step]
                'mapper' : argument of Hadoop Streaming's -mapper; if left
                    unspecified, use IdentityMapper
                'reducer' : argument of Hadoop Streaming's -reducer; if left
                    unspecified, use IdentityReducer
                'inputs' : list of input directories
                'no_input_prefix' : key that's present iff intermediate dir
                    should not be prepended to inputs
                'output' : output directory
                'no_output_prefix' : key that's present iff intermediate dir
                    should not be prepended to output dir
                'partition' : sort-like options to pass to
                    KeyFieldBasedPartitioner; of form -k pos1[,pos2], where pos
                    is ordinarily of the form f[.c][opts], but only f[opts]
                    are permitted by Dooplicity's EMR simulator
                'sort' : sort-like options to pass to sort comparator
                    KeyFieldBasedComparator; of any form permitted by UNIX
                    sort. (Hadoop only recognizes n and r though.)
                'tasks' : exact number of reduce tasks OR 'NUMx' to obtain
                    NUM * number of "slots" tasks (typically number of
                        processing cores) OR NUM1,NUM2 to decide
                    "optimal" number of tasks between NUM1 and NUM2 inclusive
                'inputformat' : input format; present only if necessary
                'archives' : archives parameter; present only if necessary
                'files' : files parameter; present only if necessary
                'multiple_outputs' : key that's present iff there are multiple
                    outputs
                'index_output' : key that's present iff output LZOs should be
                    indexed after step; applicable only in Hadoop modes
                'extra_args' : list of '-D' args
            }

        protosteps: array of protosteps
        action_on_failure: action on failure to take
        jar: path to Hadoop Streaming jar
        step_dir: where to find Python scripts for steps
        reducer_count: number of reducers; determines number of tasks
        unix: performs UNIX-like path joins; also inserts pypy in for
            executable since unix=True only on EMR
        no_direct_copy: True iff output of a step should not be written
            directly to S3; this may be turned on if there are missing-
            intermediate-data issues when speculative execution is
            also turned on

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
        identity_mapper = ('cut -f 2-' if unix else 'cat')
        final_output = (path_join(unix, intermediate_dir,
                                        protostep['output'])
                        if ('no_output_prefix' not in
                            protostep or not protostep['no_output_prefix'])
                        else protostep['output'])
        final_output_url = ab.Url(final_output)
        no_direct_copy_for_step = (
                    unix and final_output_url.is_s3 and no_direct_copy
                )
        if no_direct_copy_for_step:
            intermediate_output = _hdfs_temp_dir + final_output_url.suffix[1:]
        else:
            intermediate_output = final_output
        if 'reducer' in protostep:
            assert 'tasks' in protostep
            tasks = str(protostep['tasks'])
            assert (tasks.endswith('x') or len(tasks.split(',')) == 2 or 
                    float(tasks).is_integer())
            if tasks.endswith('x'):
                reducer_task_count = int(tasks[:-1]) * reducer_count
            else:
                tasks = tasks.split(',')
                if len(tasks) == 2:
                    min_tasks = int(tasks[0])
                    try:
                        max_tasks = int(tasks[1])
                    except ValueError:
                        max_tasks = None
                    '''Task count logic: number of reduce tasks is equal to
                    min(minimum number of tasks greater than
                        min_tasks that's a multiple of the number of reducers,
                        max_tasks); max_tasks need not be specified'''
                    if not (min_tasks % reducer_count):
                        reducer_task_count = min_tasks
                    else:
                        reducer_task_count = \
                            min_tasks + reducer_count - (
                                (min_tasks + reducer_count) % reducer_count
                            )
                    if max_tasks is not None:
                        reducer_task_count = min(reducer_task_count, max_tasks)
                else:
                    # len tasks is 1
                    reducer_task_count = int(tasks[0])
        else:
            reducer_task_count = 0
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
                                        protostep['mapper'])])
                        if 'mapper' in protostep else identity_mapper,
                reducer=' '.join(['pypy' if unix
                        else _executable, 
                        path_join(unix, step_dir,
                                        protostep['reducer'])]) 
                        if 'reducer' in protostep else 'cat',
                action_on_failure=action_on_failure,
                jar=jar,
                tasks=reducer_task_count,
                partition_options=(protostep['partition']
                    if 'partition' in protostep else None),
                sort_options=(protostep['sort']
                    if 'sort' in protostep else None),
                archives=(protostep['archives']
                    if 'archives' in protostep else None),
                files=(protostep['files']
                    if 'files' in protostep else None),
                multiple_outputs=(True if 'multiple_outputs'
                        in protostep 
                        and protostep['multiple_outputs']
                        else False
                    ),
                mod_partitioner=(True if 'mod_partitioner'
                        in protostep
                        and protostep['mod_partitioner']
                        else False
                    ),
                inputformat=(protostep['inputformat']
                    if 'inputformat' in protostep else None),
                extra_args=([extra_arg.format(task_count=reducer_count)
                    for extra_arg in protostep['extra_args']]
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
        if no_direct_copy_for_step:
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
    def __init__(self, manifest, output_dir, isofrag_idx=None,
            intermediate_dir='./intermediate', force=False, aws_exe=None,
            profile='default', region=None, service_role=None,
            instance_profile=None, verbose=False, curl_exe=None,
            max_task_attempts=4, dbgap_key=None
        ):
        '''Store all errors uncovered in a list, then output. This prevents the
        user from having to rerun Rail-RNA to find what else is wrong with
        the command-line parameters.'''
        self.errors = []
        self.manifest_dir = None
        self.manifest = manifest
        self.isofrag_idx = isofrag_idx
        self.output_dir = output_dir
        self.intermediate_dir = intermediate_dir
        self.aws_exe = aws_exe
        self.region = region
        self.force = force
        self.checked_programs = set()
        self.curl_exe = curl_exe
        self.verbose = verbose
        self.profile = profile
        self.service_role = service_role
        self.instance_profile = instance_profile
        self.dbgap_key = dbgap_key
        if not (float(max_task_attempts).is_integer()
                        and max_task_attempts >= 1):
            self.errors.append('Max task attempts (--max-task-attempts) '
                               'must be an integer greater than 0, but '
                               '{0} was entered'.format(
                                                    max_task_attempts
                                                ))
        self.max_task_attempts = max_task_attempts

    def check_s3(self, reason=None, is_exe=None, which=None):
        """ Checks for AWS CLI and configuration file.

            In this script, checking is performed as soon as it is found
            that thE CLI is needed. If anything is awry, a RuntimeError is
            raised _immediately_ (the standard behavior is to raise a
            RuntimeError only after errors are accumulated). A reason
            specifying where credentials were first needed can also be
            provided.

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
                        if tokens[0] == 'region':
                            if self.region is None:
                                self.region = tokens[1]
                            region_holder = tokens[1]
                        elif tokens[0] == 'aws_access_key_id':
                            self._aws_access_key_id = tokens[1]
                        elif tokens[0] == 'aws_secret_access_key':
                            self._aws_secret_access_key = tokens[1]
                        elif tokens[0] == 'emr':
                            grab_roles = True
                        elif tokens[0] == 'service_role' and grab_roles:
                            self.service_role = tokens[1]
                        elif tokens[0] == 'instance_profile' and grab_roles:
                            self.instance_profile = tokens[1]
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
                                    'To set this file up, run "aws configure" '
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
        if self.region is None:
            self.region = 'us-east-1'
        # Finalize roles; if they're still None, try default Amazonian ones
        if self.service_role is None:
            print_to_screen('Warning: IAM service role not found. Attempting '
                            '"EMR_DefaultRole".',
                            newline=True, carriage_return=False)
            self.service_role = 'EMR_DefaultRole'
        if self.instance_profile is None:
            print_to_screen('Warning: EC2 instance profile not found. '
                            'Attempting "EMR_EC2_DefaultRole".',
                            newline=True, carriage_return=False)
            self.instance_profile = 'EMR_EC2_DefaultRole'
        self.checked_programs.add('AWS CLI')

    def check_cloudformation(self, stack_name):
        """ Gets subnet ID, security groups, roles, etc. from secure stack.

            This is basically a version of check_s3() above tailored
            for describing a stack and setting appropriate object members.

            stack_name: name of stack
            region: region in which stack was created

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
                        if tokens[0] == 'region':
                            if self.region is None:
                                self.region = tokens[1]
                            region_holder = tokens[1]
                        elif tokens[0] == 'aws_access_key_id':
                            self._aws_access_key_id = tokens[1]
                        elif tokens[0] == 'aws_secret_access_key':
                            self._aws_secret_access_key = tokens[1]
                        elif tokens[0] == 'emr':
                            grab_roles = True
                        elif tokens[0] == 'service_role' and grab_roles:
                            self.service_role = tokens[1]
                        elif tokens[0] == 'instance_profile' and grab_roles:
                            self.instance_profile = tokens[1]
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
                                    'To set this file up, run "aws configure" '
                                    'after installing the AWS CLI.')
                                )
        if len(self.errors) != original_errors_size:
            raise RuntimeError(('\n'.join(['%d) %s' % (i+1, error)
                                for i, error
                                in enumerate(self.errors)]) 
                                if len(self.errors) > 1
                                else self.errors[0]) + 
                                '\n\nNote that the AWS CLI is needed '
                                'because you are attempting to launch '
                                'Rail-RNA into a VPC from a secure '
                                'stack.')
        if self.region is None:
            self.region = 'us-east-1'
        # Get appropriate vars
        command_to_run = ('{} cloudformation describe-stacks '
                                     '--stack-name {} --region {}').format(
                                            self.aws_exe, stack_name, region
                                        )
        try:
            described_stack = subprocess.check_output(command_to_run,
                                                        shell=True,
                                                        executable='/bin/bash')
        except subprocess.CalledProcessError:
            self.errors.append('Error encountered attempting to run AWS CLI '
                               'to describe stack {}.'.format(stack_name))
            raise_runtime_error(self)
        else:
            if 'does not exist' in described_stack:
                base.errors.append((
                        'Stack name "{}" does not exist in region "{}".'
                    ).format(stack_name, region))
                raise_runtime_error(self)
            elif 'CREATE_COMPLETE' not in described_stack:
                self.errors.append((
                        'State of stack "{stack_name}" is not '
                        'CREATE_COMPLETE. Stack state '
                        'must be CREATE_COMPLETE before proceeding; if you '
                        'just initiated stack creation with the AWS CLI, you '
                        'can check whether it\'s done using '
                        '"aws cloudformation describe-stacks --stack-name '
                        '{stack_name} --region {region}".'
                    ).format(stack_name=stack_name, region=self.region))
                raise_runtime_error(self)
            else:
                described_stack = json.loads(
                                        described_stack
                                    )['Stacks'][0]['Outputs']
                outputs = {}
                for output in described_stack:
                    outputs[output['OutputKey']] = output['OutputValue']
                for output in ['ServiceRole', 'InstanceProfile',
                                'MasterSecurityGroupId',
                                'SlaveSecurityGroupId', 'PublicSubnetId',
                                'SecureBucketName']
                    if output not in outputs:
                        self.errors.append(('{} not among outputs when '
                                           'describing stack {}. '
                                           'Make sure you used an approved '
                                           'CloudFormation template to create '
                                           'the stack.').format(output,
                                                                stack_name))
                raise_runtime_error(self)
                self.service_role = outputs['ServiceRole']
                self.instance_profile = outputs['InstanceProfile']
                self.ec2_subnet_id = outputs['PublicSubnetId']
                self.ec2_master_security_group_id = outputs[
                                                        'MasterSecurityGroupId'
                                                    ]
                self.ec2_slave_security_group_id = outputs[
                                                        'SlaveSecurityGroupId'
                                                    ]
                self.secure_bucket_name = outputs['SecureBucketName']
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
            reason: FOR CURL ONLY: raise RuntimeError _immediately_ if cURL
                not found but needed
            is_exe: is_exe function
            which: which function

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
            which_entered_exe = which(entered_exe)
            if which_entered_exe is None:
                self.errors.append(
                    ('The executable "{0}" entered for {1} via {2} was '
                     'either not found or is not executable.').format(exe,
                                                                program_name,
                                                                parameter)
                )
            to_return = which_entered_exe
        else:
            to_return = entered_exe
        if original_errors_size != len(self.errors) and reason:
            raise RuntimeError((('\n'.join(['%d) %s' % (i+1, error)
                                for i, error
                                in enumerate(self.errors)])
                                if len(self.errors) > 1 else self.errors[0]) + 
                                '\n\nNote that cURL is needed because {0}.'
                                ' If all dependence on web resources is '
                                'removed from the pipeline, cURL need '
                                'not be installed.').format(reason))
        self.checked_programs.add(program_name)
        return to_return

    @staticmethod
    def add_args(general_parser, exec_parser, required_parser):
        exec_parser.add_argument(
            '--aws', type=str, required=False, metavar='<exe>',
            default=exe_paths.aws,
            help=('path to AWS CLI executable (def: %s)'
                    % (exe_paths.aws
                        if exe_paths.aws is not None
                        else 'aws'))
        )
        exec_parser.add_argument(
            '--curl', type=str, required=False, metavar='<exe>',
            default=exe_paths.curl,
            help=('path to cURL executable (def: %s)'
                    % (exe_paths.curl
                        if exe_paths.curl is not None
                        else 'curl'))
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
            help='Myrna-style manifest file'
        )
        '''--region and --max-task-attempts help looks different from mode to '
        ' mode; don't include them here. Also postpone inclusion of isofrag
        index; should be an in algo parser.'''

def raise_runtime_error(bases):
    """ Raises RuntimeError if any base.errors is nonempty.

        bases: dictionary mapping IPython engine IDs to RailRnaErrors
            instances or a single RailRnaError instance.

        No return value.
    """
    assert isinstance(bases, RailRnaErrors) or isinstance(bases, dict)
    if isinstance(bases, RailRnaErrors) and bases.errors:
        raise RuntimeError(
                '\n'.join(
                        ['%d) %s' % (i+1, error) for i, error
                            in enumerate(bases.errors)]
                    ) if len(bases.errors) > 1 else bases.errors[0]
            )
    elif isinstance(bases, dict):
        errors_to_report = defaultdict(set)
        for engine in bases:
            assert isinstance(bases[engine], RailRnaErrors)
            if bases[engine].errors:
                errors_to_report['\n'.join(
                            ['%d) %s' % (i+1, error) for i, error
                                in enumerate(bases[engine].errors)]
                        ) if len(bases[engine].errors) > 1
                          else bases[engine].errors[0]].add(engine)
        runtimeerror_message = []
        if errors_to_report:
            for message in errors_to_report:
                runtimeerror_message.extend(
                    ['Engine(s) %s report(s) the following errors.'
                        % engine_string_from_list(
                              errors_to_report[message]
                            ), message]
                    )
            raise RuntimeError('\n',join(errors_to_report))

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
            rc = Client(profile=ipython_profile)
        except ValueError:
            errors.append(
                    'Cluster configuration profile "%s" was not '
                    'found.' % ipython_profile
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
        except UnboundLocalError:
            # Client referenced before assignment; arises from ImportError
            pass
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
    print_to_screen('Detected %d running IPython engines.' 
                                % len(rc.ids))
    # Use Dill to permit general serializing
    try:
        import dill
    except ImportError:
        raise RuntimeError(
                'Rail-RNA requires Dill in "parallel" mode. Install it by '
                'running "pip install dill", or see the StackOverflow '
                'question http://stackoverflow.com/questions/23576969/'
                'how-to-install-dill-in-ipython for other leads.'
            )
    else:
        rc[:].use_dill()
    return rc

class RailRnaLocal(object):
    """ Checks local- or parallel-mode JSON for programs and input parameters.

        Subsumes only those parameters relevant to local mode. Adds errors
        to base instance of RailRnaErrors.
    """
    def __init__(self, base, check_manifest=False,
                    num_processes=1, keep_intermediates=False,
                    gzip_intermediates=False, gzip_level=3,
                    sort_memory_cap=(300*1024), parallel=False,
                    local=True, scratch=None, direct_write=False,
                    ansible=None, do_not_copy_index_to_nodes=False,
                    sort_exe=None, fastq_dump_exe=None, vdb_config_exe=None):
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
        base.do_not_copy_index_to_nodes = do_not_copy_index_to_nodes
        base.direct_write = direct_write
        if not parallel:
            if output_dir_url.is_local:
                if os.path.exists(output_dir_url.to_url()):
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
            # Check transcript fragment index
            if base.isofrag_idx is not None:
                isofrag_url = ab.Url(base.isofrag_idx)
                if (isofrag_url.is_s3 and 'AWS CLI'
                        not in base.checked_programs):
                    base.check_s3(reason='the isofrag index is on S3',
                                    is_exe=is_exe,
                                    which=which)
                    # Change ansible params
                    ansible.aws_exe = base.aws_exe
                    ansible.profile = base.profile
                elif isofrag_url.is_curlable \
                    and 'cURL' not in base.checked_programs:
                    base.curl_exe = base.check_program(
                                    'curl', 'cURL', '--curl',
                                    entered_exe=base.curl_exe,
                                    reason='the isofrag index is on the web',
                                    is_exe=is_exe,
                                    which=which
                                )
                    ansible.curl_exe = base.curl_exe
                if not ansible.exists(isofrag_url.to_url()):
                    base.errors.append(('Isofrag index (--isofrag-idx) {0} '
                                        'does not exist. Check the URL and '
                                        'try again.').format(base.isofrag_idx))
                else:
                    if not isofrag_url.is_local:
                        '''Download isofrag index only if not an IPython engine
                        (not parallel).'''
                        base.isofrag_dir = base.intermediate_dir
                        base.isofrag_idx = os.path.join(
                                            base.isofrag_dir,
                                            os.path.basename(base.isofrag_idx)
                                        )
                        ansible.get(isofrag_url.to_url(),
                                    destination=base.isofrag_idx)
                    base.isofrag_idx = os.path.abspath(base.isofrag_idx)
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
                and 'cURL' not in base.checked_programs:
                base.curl_exe = base.check_program('curl', 'cURL', '--curl',
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
                    base.manifest = os.path.join(base.manifest_dir, 'MANIFEST')
                    ansible.get(manifest_url.to_url(),
                                destination=base.manifest)
                base.manifest = os.path.abspath(base.manifest)
                files_to_check = []
                base.sample_count = 0
                if base.dbgap_key is not None \
                    and not os.path.exists(dbgap_key):
                    base.errors.append(('dbGaP repository key file '
                                        '(--dbgap-key) "{}" '
                                        'does not exist. Check its path and '
                                        'try again. Note that it must be '
                                        'stored on the local '
                                        'filesystem.').format(dbgap_key))
                base.dbgap_present = False
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
                            single_url = ab.Url(tokens[0])
                            if single_url.is_dbgap:
                                base.dbgap_present = True
                        else:
                            base.errors.append(('The following line from the '
                                                'manifest file {0} '
                                                'has an invalid number of '
                                                'tokens:\n{1}'
                                                ).format(
                                                        manifest_url.to_url(),
                                                        line.strip()
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
                                                        line.strip()
                                                    ))
                if base.dbgap_present:
                    base.errors.append('Rail-RNA does not currently work '
                                       'with dbGaP accession numbers '
                                       'in manifest files in local and '
                                       'parallel modes. Decrypt your data '
                                       'first and try again with FASTQ paths '
                                       'in the manifest file instead.')
                '''if base.dbgap_present and base.dbgap_key is None:
                    base.errors.append('dbGaP accession numbers are in '
                                       'manifest file, but no dbGaP '
                                       'repository key file (--dbgap-key) was '
                                       'provided. ')'''
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
                                and 'cURL' not in base.checked_programs:
                                if local:
                                    base.curl_exe = base.check_program('curl',
                                                    'cURL',
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
                            elif filename_url.is_sra:
                                if 'SRA Tools' not in base.checked_programs:
                                    base.fastq_dump_exe = base.check_program(
                                                'fastq-dump', 'fastq-dump',
                                                '--fastq-dump',
                                                entered_exe=fastq_dump_exe,
                                                is_exe=is_exe,
                                                which=which
                                            )
                                    base.vdb_config_exe = base.check_program(
                                                'vdb-config', 'vdb-config',
                                                '--vdb-config',
                                                entered_exe=vdb_config_exe,
                                                is_exe=is_exe,
                                                which=which
                                            )
                                    fastq_dump_version_command = [
                                            base.fastq_dump_exe, '--version'
                                        ]
                                    vdb_config_version_command = [
                                            base.vdb_config_exe, '--version'
                                        ]
                                    try:
                                        fastq_dump_version = \
                                            subprocess.check_output(
                                                fastq_dump_version_command
                                            ).strip().split(':')[1].strip()
                                    except Exception as e:
                                        base.errors.append(
                                            ('Error "{0}" encountered '
                                             'attempting to execute '
                                             '"{1}".').format(
                                                    e.message,
                                                    ' '.join(
                                                    fastq_dump_version_command
                                                   )
                                                ))
                                    else:
                                        try:
                                            vdb_config_version = \
                                                subprocess.check_output(
                                                    vdb_config_version_command
                                                ).strip().split(':')[1].strip()
                                        except Exception as e:
                                            base.errors.append(
                                                ('Error "{0}" encountered '
                                                 'attempting to execute '
                                                 '"{1}".').format(
                                                    e.message,
                                                    ' '.join(
                                                    vdb_config_version_command
                                                   )
                                                ))
                                        else:
                                            if fastq_dump_version \
                                                != vdb_config_version:
                                                base.errors.append(
                                                    ('fastq-dump and '
                                                     'vdb-config '
                                                     'versions do not agree: '
                                                     'fastq-dump '
                                                     'version is {}, while '
                                                     'vdb-config version '
                                                     'is {}. '
                                                     'Reinstall '
                                                     'SRA Tools v2.5 or '
                                                     'greater and '
                                                     'try again.').format(
                                                            fastq_dump_version,
                                                            vdb_config_version
                                                        )
                                                    )
                                            base.sra_tools_version \
                                                = fastq_dump_version
                                    base.checked_programs.add('SRA Tools')
                            if not filename_url.is_sra \
                                and not ansible.exists(filename):
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
            raise_runtime_error(base)
            from multiprocessing import cpu_count
            if num_processes:
                if not (float(num_processes).is_integer()
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
                if not (float(gzip_level).is_integer()
                                        and 9 >= gzip_level >= 1):
                    base.errors.append('Gzip level (--gzip-level) '
                                       'must be an integer between 1 and 9, '
                                       'but {0} was entered.'.format(
                                                        gzip_level
                                                    ))
            base.gzip_intermediates = gzip_intermediates
            base.gzip_level = gzip_level
            if not (sort_memory_cap > 0):
                base.errors.append('Sort memory cap (--sort-memory-cap) '
                                   'must take a value larger than 0, '
                                   'but {0} was entered.'.format(
                                                        sort_memory_cap
                                                    ))
            base.sort_memory_cap = sort_memory_cap
        if parallel and scratch:
            expanded_scratch = os.path.expandvars(scratch)
            if not os.path.exists(expanded_scratch):
                try:
                    os.makedirs(expanded_scratch)
                except OSError:
                    base.errors.append(
                            ('Could not create scratch directory %s; '
                             'check that it\'s not a file and that '
                             'write permissions are active.')
                                % expanded_scratch
                        )
        base.scratch = scratch
        if sort_exe:
            sort_exe_parameters = [parameter.strip()
                                    for parameter in sort_exe.split(' ')]
        else:
            sort_exe_parameters = []
        check_scratch = True
        try:
            sort_scratch = sort_exe_parameters[
                    sort_exe_parameters.index('--temporary-directory')+1
                ]
        except IndexError:
            base.errors.append(
                    ('"--temporary-directory" parameter was passed '
                     'to sort executable without specifying temporary '
                     'directory')
                )
        except ValueError:
            try:
                sort_scratch = sort_exe_parameters[
                        sort_exe_parameters.index('-T')+1
                    ]
            except IndexError:
                base.errors.append(
                    ('"-T" parameter was passed '
                     'to sort executable without specifying temporary '
                     'directory')
                )
            except ValueError:
                sort_scratch = base.scratch
                check_scratch = False
        if parallel and check_scratch:
            expanded_sort_scratch  = os.path.expandvars(sort_scratch)
            if not os.path.exists(expanded_sort_scratch):
                try:
                    os.makedirs(expanded_sort_scratch)
                except OSError:
                    base.errors.append(
                            ('Could not create sort scratch directory %s; '
                             'check that it\'s not a file and that '
                             'write permissions are active.')
                            % expanded_sort_scratch
                        )
        base.sort_exe = ' '.join(
                            [base.check_program('sort', 'sort', '--sort',
                                entered_exe=(sort_exe_parameters[0]
                                                if sort_exe_parameters
                                                else None),
                                is_exe=is_exe,
                                which=which)] + sort_exe_parameters[1:]
                             + (['-T', base.scratch] if (not check_scratch
                                 and base.scratch is not None) else [])
                        )

    @staticmethod
    def add_args(required_parser, general_parser, output_parser, 
                    exec_parser, prep=False, align=False, parallel=False):
        """ Adds parameter descriptions relevant to local mode to an object
            of class argparse.ArgumentParser.

            prep: preprocess-only
            align: align-only
            parallel: add parallel-mode arguments

            No return value.
        """
        exec_parser.add_argument(
            '--sort', type=str, required=False, metavar='<exe>',
            default=exe_paths.sort,
            help=('path to sort executable; include extra sort parameters '
                  'here (def: %s)'
                    % (exe_paths.sort
                        if exe_paths.sort is not None else 'sort'))
        )
        if align:
            required_parser.add_argument(
                '-i', '--input', type=str, required=True, metavar='<dir>',
                help='input directory with preprocessed reads; must be local'
            )
        '''else:
            # "prep" or "go" flows
            general_parser.add_argument(
                '--dbgap-key', required=False, metavar='<file>'
                default=None,
                help='path to dbGaP key file, which has the extension "ngc"; '
                     'must be on local filesystem'
            )'''
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
            general_parser.add_argument(
                '--scratch', type=str, required=False, metavar='<dir>',
                default=None,
                help=('directory for storing temporary files; BASH variables '
                      'specified with dollar signs are recognized here '
                      '(def: securely created temporary directory)')
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
                      'profile)')
            )
            general_parser.add_argument(
                '--scratch', type=str, required=False, metavar='<dir>',
                default=None,
                help=('scratch directory for storing '
                      'Bowtie index and temporary files before they are '
                      'committed; node- or slot-local BASH variables '
                      'specified with escaped dollar signs are recognized '
                      'here (def: return value of tempfile.gettempdir(); see '
                      'Python docs)')
            )
            general_parser.add_argument(
                '--direct-write', action='store_const', const=True,
                default=False,
                help=('write intermediate files directly to the log directory '
                      'rather than first writing to scratch and moving '
                      'results')
            )
            if not prep:
                general_parser.add_argument(
                        '--do-not-copy-index-to-nodes', action='store_const',
                        const=True,
                        default=False,
                        help=('does not copy Bowtie/Bowtie 2 indexes to '
                              'nodes before starting job flow; copying '
                              'requires Herd')
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
        general_parser.add_argument(
            '-r', '--sort-memory-cap', type=float, required=False,
            metavar='<dec>',
            default=(300*1024),
            help=('maximum amount of memory (in bytes) used by UNIX sort '
                  'per process')
        )
        general_parser.add_argument(
            '--max-task-attempts', type=int, required=False,
            metavar='<int>',
            default=(4 if parallel else 1),
            help=('maximum number of attempts per task')
        )

class RailRnaElastic(object):
    """ Checks elastic-mode input parameters and relevant programs.

        Subsumes only those parameters relevant to elastic mode. Adds errors
        to base instance of RailRnaErrors.
    """
    def __init__(self, base, check_manifest=False,
        log_uri=None, ami_version='3.8.0',
        visible_to_all_users=False, tags='',
        name='Rail-RNA Job Flow',
        action_on_failure='TERMINATE_JOB_FLOW',
        hadoop_jar=None,
        master_instance_count=1, master_instance_type='c1.xlarge',
        master_instance_bid_price=None, core_instance_count=1,
        core_instance_type=None, core_instance_bid_price=None,
        task_instance_count=0, task_instance_type=None,
        task_instance_bid_price=None, ec2_key_name=None,
        ec2_subnet_id=None, ec2_master_security_group_id=None,
        ec2_slave_security_group_id=None, keep_alive=False,
        termination_protected=False, consistent_view=False,
        no_direct_copy=False, intermediate_lifetime=4,
        secure_stack_name=None):

        base.secure_stack_name = secure_stack_name
        # CLI is REQUIRED in elastic mode
        if base.secure_stack_name is None:
            base.check_s3(reason='Rail-RNA is running in "elastic" mode')
        else:
            # CloudFormation version required
            base.check_cloudformation(stack_name=base.secure_stack_name)

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
            if not (float(intermediate_lifetime).is_integer() and
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
        # Set directory for storing Rail-RNA, bootstraps, and possibly manifest
        base.dependency_dir = base.output_dir + '.dependencies'
        dependency_dir_url = ab.Url(base.dependency_dir)
        if dependency_dir_url.is_s3:
            # Set up rule on S3 for deleting dependency dir
            final_dependency_dir = dependency_dir_url.to_url() + '/'
            while final_dependency_dir[-2] == '/':
                final_dependency_dir = final_dependency_dir[:-1]
            ansible.s3_ansible.expire_prefix(final_dependency_dir,
                                                days=intermediate_lifetime)
        # Check manifest; download it if necessary
        manifest_url = ab.Url(base.manifest)
        if manifest_url.is_curlable \
            and 'cURL' not in base.checked_programs:
            base.curl_exe = base.check_program('curl', 'cURL', '--curl',
                                    entered_exe=base.curl_exe,
                                    reason='the manifest file is on the web')
            ansible.curl_exe = base.curl_exe
        if not ansible.exists(manifest_url.to_url()):
            base.errors.append(('Manifest file (--manifest) {0} '
                                'does not exist. Check the URL and '
                                'try again.').format(base.manifest))
        if base.isofrag_idx is not None:
            isofrag_url = ab.Url(base.isofrag_idx)
            if isofrag_url.is_curlable \
                and 'cURL' not in base.checked_programs:
                base.curl_exe = base.check_program(
                                    'curl', 'cURL', '--curl',
                                    entered_exe=base.curl_exe,
                                    reason='the isofrag index is on the web'
                                )
                ansible.curl_exe = base.curl_exe
            if not ansible.exists(isofrag_url.to_url()):
                base.errors.append(('Isofrag index (--isofrag-idx) {0} '
                                    'does not exist. Check the URL and '
                                    'try again.').format(base.isofrag_idx))
        raise_runtime_error(base)
        if not manifest_url.is_local:
            temp_manifest_dir = tempfile.mkdtemp()
            from tempdel import remove_temporary_directories
            register_cleanup(remove_temporary_directories,
                                [temp_manifest_dir])
            manifest = os.path.join(temp_manifest_dir, 'MANIFEST')
            ansible.get(base.manifest, destination=manifest)
        else:
            manifest = manifest_url.to_url()
        if base.dbgap_key is not None and \
            not os.path.exists(base.dbgap_key):
            base.errors.append(('dbGaP repository key file '
                                '(--dbgap-key) "{}" '
                                'does not exist. Check its path and '
                                'try again. Note that it must be '
                                'found on the local '
                                'filesystem; Rail-RNA uploads this '
                                'file to S3 securely with server-side '
                                'encryption enabled and schedules '
                                'it for deletion.').format(dbgap_key))
        base.dbgap_present = False
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
                    single_url = ab.Url(tokens[0])
                    if single_url.is_dbgap:
                        base.dbgap_present = True
                        sra_tools_needed = True
                    elif single_url.is_sra:
                        sra_tools_needed = True
                else:
                    check_sample_label = False
                    base.errors.append(('The following line from the '
                                        'manifest file {0} '
                                        'has an invalid number of '
                                        'tokens:\n{1}'
                                        ).format(
                                                manifest_url.to_url(),
                                                line.strip()
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
                                                line.strip()
                                            ))
        if base.dbgap_present and base.dbgap_key is None:
            base.errors.append('dbGaP accession numbers are in '
                               'manifest file, but no dbGaP '
                               'repository key file (--dbgap-key) was '
                               'provided. ')
        if base.dbgap_present and base.secure_stack_name is None:
            # Secure stack must be specified!
            base.errors.append('Manifest file (--manifest) has dbGaP data, '
                               'but secure stack name (--secure-stack-name) '
                               'was not specified. Rail-RNA must be launched '
                               'into a VPC associated with some secure stack '
                               'in elastic mode when analyzing dbGaP data.')
            raise_runtime_error(base)
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
                        and 'cURL' not in base.checked_programs:
                        base.curl_exe = base.check_program('curl', 'cURL',
                                            '--curl',
                                            entered_exe=base.curl_exe,
                                            reason=('at least one sample '
                                              'FASTA/FASTQ from the '
                                              'manifest file is on '
                                              'the web'))
                        ansible.curl_exe = base.curl_exe
                    if not filename_url.is_sra \
                        and not ansible.exists(filename_url.to_url()):
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
        base.ec2_subnet_id = ec2_subnet_id
        base.ec2_master_security_group_id = ec2_master_security_group_id
        base.ec2_slave_security_group_id = ec2_slave_security_group_id
        # Bail before copying anything to S3
        raise_runtime_error(base)
        if base.secure_stack_name is not None:
            if (ec2_subnet_id is None
                    or base.ec2_master_security_group_id is None
                    or base.ec2_slave_security_group_id is None):
                raise RuntimeError('Manifest file (--manifest) includes dbGaP '
                                   'data and/or Rail-RNA is being run in '
                                   'secure mode, so the EMR cluster must be '
                                   'launched into a public VPC subnet '
                                   'specified with --ec2-subnet-id, where '
                                   'appropriate security groups are specified '
                                   'with --ec2-master-security-group-id and '
                                   '--ec2-slave-security-group-id . Google '
                                   'for "NIH Best Practices for '
                                   'Controlled-Access Data" '
                                   'and the Whalley-Pizarro whitepaper '
                                   '"Architecting for Genomic Data '
                                   'Security and Compliance in AWS" for '
                                   'more information. CloudFormation '
                                   'templates in '
                                   '$RAILDOTBIO/rail-rna/cloudformation '
                                   'can create the required stacks. Refer '
                                   'to the Rail documentation for '
                                   'instructions on their proper use.')
            else:
                if (ansibles.bucket_from_url(base.intermediate_dir)
                    != base.secure_bucket_name):
                    print_to_screen(('Warning: intermediate directory '
                                       '"{}" not in secure bucket "{}", which '
                                       'is associated with secure stack '
                                       '"{}".').format(base.intermediate_dir,
                                                    base.secure_bucket_name,
                                                    base.secure_stack_name),
                                        newline=True, carriage_return=True,
                                    )
                if (ansibles.bucket_from_url(base.output_dir)
                    != base.secure_bucket_name):
                    print_to_screen(('Warning: output directory '
                                       '"{}" not in secure bucket "{}", which '
                                       'is associated with secure stack '
                                       '"{}".').format(base.output_dir,
                                                    base.secure_bucket_name,
                                                    base.secure_stack_name),
                                        newline=True, carriage_return=True,
                                    )
                question = ('Manifest file (--manifest) includes dbGaP data '
                            'and/or Rail-RNA is being run in secure mode. '
                            'Do you certify that the EC2 subnet ID '
                            '"{ec2_subnet_id}", EC2 master security group ID '
                            '"{ec2_master_security_group_id}", and EC2 slave '
                            'security group ID '
                            '"{ec2_slave_security_group_id}" '
                            'correspond to a stack that provides '
                            'the security features included in one of the '
                            'dbGaP CloudFormation templates in '
                            '$RAILDOTBIO/rail-rna/cloudformation?').format(
                                    ec2_subnet_id=base.ec2_subnet_id,
                                    ec2_master_security_group_id=\
                                        base.ec2_master_security_group_id,
                                    ec2_slave_security_group_id=\
                                        base.ec2_slave_security_group_id
                                )
                while True:
                    sys.stdout.write('%s [y/n]: ' % question)
                    try:
                        try:
                            if not strtobool(raw_input().lower()):
                                print_to_screen(
                                    'Set up a secure VPC and try again.',
                                    newline=True, carriage_return=True,
                                )
                            else: break
                        except KeyboardInterrupt:
                            sys.stdout.write('\n')
                            sys.exit(0)
                    except ValueError:
                        sys.stdout.write('Please enter \'y\' or \'n\'.\n')
        # Download+upload isofrag index if necessary
        if base.isofrag_idx is not None and isofrag_url.is_curlable:
            temp_isofrag_dir = tempfile.mkdtemp()
            from tempdel import remove_temporary_directories
            register_cleanup(remove_temporary_directories,
                                [temp_isofrag_dir])
            isofrag = os.path.join(temp_isofrag_dir,
                                    os.path.basename(base.isofrag_idx))
            ansible.get(base.isofrag_idx, destination=isofrag)
            base.isofrag_idx = path_join(
                                    True,
                                    base.output_dir + '.isofrag',
                                    os.path.basename(base.isofrag_idx)
                                )
            ansible.put(isofrag, base.isofrag_idx)
            shutil.rmtree(temp_isofrag_dir, ignore_errors=True)

        # Upload NGC file to S3
        if base.dbgap_present:
            base.dbgap_s3_path = path_join(
                                    True,
                                    base.intermediate_dir,
                                    os.path.basename(base.dbgap_key)
                                )
            ansible.put(base.dbgap_key, base.dbgap_s3_path)

        if not manifest_url.is_s3 and output_dir_url.is_s3:
            # Copy manifest file to S3 before job flow starts
            base.manifest = path_join(True, base.dependency_dir,
                                            'MANIFEST')
            ansible.put(manifest, base.manifest)
        if not manifest_url.is_local:
            # Clean up
            shutil.rmtree(temp_manifest_dir)
        print_to_screen('Copying Rail-RNA and bootstraps to S3...',
                        newline=False, carriage_return=True)
        # Create and copy bootstraps to S3; compress and copy Rail-RNA to S3
        temp_dependency_dir = tempfile.mkdtemp()
        from tempdel import remove_temporary_directories
        register_cleanup(remove_temporary_directories, [temp_dependency_dir])
        copy_index_bootstrap = os.path.join(temp_dependency_dir,
                                                'install-index.sh')
        with open(copy_index_bootstrap, 'w') as script_stream:
             print >>script_stream, (
"""#!/usr/bin/env bash
set -e

mkdir -p $1
cd $1

fn=`basename $2`
echo fn
echo $1
echo $(pwd)
cat >.download_idx.py <<EOF
\"""
.download_idx.py

Checks if file size hasn't changed while downloading idx, and if it hasn't,
kills subprocess and restarts it.
\"""

import sys
import time
import os
import subprocess

tries = 0
url = sys.argv[1]
filename = sys.argv[2]
while tries < 5:
    break_outer_loop = False
    s3cmd_process \\
        = subprocess.Popen(['s3cmd', 'get', url, './', '-f',
                    '--add-header="x-amz-request-payer: requester"'])
    time.sleep(1)
    last_check_time = time.time()
    try:
        last_size = os.path.getsize(filename)
    except OSError:
        last_size = 0
    while s3cmd_process.poll() is None:
        now_time = time.time()
        if now_time - last_check_time > 120:
            try:
                new_size = os.path.getsize(filename)
            except OSError:
                new_size = 0
            if new_size == last_size:
                # Download stalled
                break_outer_loop = True
                break
            else:
                last_size = new_size
            last_check_time = now_time
            time.sleep(1)
    if break_outer_loop:
        tries += 1
        s3cmd_process.kill()
        try:
            os.remove(filename)
        except OSError:
            pass
        time.sleep(2)
        continue
    if s3cmd_process.poll() > 0:
        if tries > 5:
            raise RuntimeError(('Non-zero exitlevel %d from s3cmd '
                                'get command')  % (
                                                s3cmd_process.poll()
                                            ))
        else:
            try:
                os.remove(filename)
            except OSError:
                pass
            tries += 1
            time.sleep(2)
            continue                       
    break
if tries > 5:
    raise RuntimeError('Could not download file from S3 '
                       'within 5 tries.')

EOF

python .download_idx.py $2 $fn
tar xvzf $fn
cd index
ln -s genome.4.bt2 genome.4.ebwt
cd ..
"""
                )
        base.copy_index_bootstrap = path_join(
                                        True, base.dependency_dir,
                                        os.path.basename(copy_index_bootstrap)
                                    )
        ansible.put(copy_index_bootstrap, base.copy_index_bootstrap)
        os.remove(copy_index_bootstrap)
        install_s3cmd_bootstrap = os.path.join(temp_dependency_dir,
                                                'install-s3cmd.sh')
        with open(install_s3cmd_bootstrap, 'w') as script_stream:
             print >>script_stream, (
"""#!/usr/bin/env bash
set -e
export HOME=/home/hadoop
printf '\\nexport HOME=/home/hadoop\\n' >>/home/hadoop/conf/hadoop-user-env.sh
sudo ln -s /home/hadoop/.s3cfg /home/.s3cfg

curl -OL https://github.com/s3tools/s3cmd/releases/download/v1.5.2/s3cmd-1.5.2.tar.gz
tar xvzf s3cmd-1.5.2.tar.gz
cd s3cmd-1.5.2
sudo ln -s $(pwd)/s3cmd /usr/bin/s3cmd

cat >~/.s3cfg <<EOF
[default]
access_key = 
secret_key = 
security_token = 
socket_timeout = 300
multipart_chunk_size_mb = 15
reduced_redundancy = False
send_chunk = 4096
use_https = True
EOF
"""
                )
        base.install_s3cmd_bootstrap = path_join(
                                        True, base.dependency_dir,
                                        os.path.basename(
                                                install_s3cmd_bootstrap
                                            )
                                    )
        ansible.put(install_s3cmd_bootstrap, base.install_s3cmd_bootstrap)
        os.remove(install_s3cmd_bootstrap)
        # Zip Rail and put it on S3
        rail_zipped = os.path.join(temp_dependency_dir, 'rail-rna.zip')
        target_to_zip = os.path.join(temp_dependency_dir, 'src')
        shutil.copytree(base_path, target_to_zip,
                        ignore=shutil.ignore_patterns('*.pyc', '.DS_Store'))
        with cd(temp_dependency_dir):
            shutil.make_archive(rail_zipped[:-4], 'zip', 'src')
        base.elastic_rail_path = path_join(
                                        True, base.dependency_dir,
                                        os.path.basename(
                                                rail_zipped
                                            )
                                    )
        ansible.put(rail_zipped, base.elastic_rail_path)
        shutil.rmtree(target_to_zip)
        install_rail_bootstrap = os.path.join(temp_dependency_dir,
                                                'install-rail.sh')
        with open(install_rail_bootstrap, 'w') as script_stream:
             print >>script_stream, (
"""#!/usr/bin/env bash
# Installs Rail-RNA and compiles required classes
# $1: where to find Rail
# $2: where JARs should be copied
# $3 on: args to pass to Rail-RNA installer
set -e
export HOME=/home/hadoop

RAILZIP=$1
JARTARGET=$2
mkdir -p ${{JARTARGET}}
shift 2
s3cmd get $RAILZIP
mkdir -p sandbox
cd sandbox
unzip ../{rail_zipped} -d ./
cd hadoop
for JAR in relevant-elephant multiple-files mod-partitioner
do
    rm -rf ${{JAR}}_out
    mkdir -p ${{JAR}}_out
    javac -classpath $(find ~/share/ *.jar | tr '\\n' ':') -d ${{JAR}}_out ${{JAR}}/*.java
    jar -cvf ${{JAR}}.jar -C ${{JAR}}_out .
    mv ${{JAR}}.jar ${{JARTARGET}}/
done
cd ../..
rm -rf sandbox
sudo python27 {rail_zipped} $@
""".format(rail_zipped=os.path.basename(rail_zipped))
            )
        base.install_rail_bootstrap = path_join(
                                        True, base.dependency_dir,
                                        os.path.basename(
                                                install_rail_bootstrap
                                            )
                                    )
        ansible.put(install_rail_bootstrap, base.install_rail_bootstrap)
        copy_bootstrap = os.path.join(temp_dependency_dir,
                                                's3cmd_s3.sh')
        base.fastq_dump_exe = _elastic_fastq_dump_exe
        base.vdb_config_exe = _elastic_vdb_config_exe
        with open(copy_bootstrap, 'w') as script_stream:
             print >>script_stream, (
"""#!/usr/bin/env bash
# s3cmd_s3.sh
#
# Download a file from an S3 bucket to given directory.  Optionally rename it.
#
# Arguments are:
# 1. s3:// URL to copy from
# 2. Local directory to copy to
# 3. If specified, name to rename file to
set -e
mkdir -p ${2}
cd ${2}
fn=`basename ${1}`
s3cmd get ${1} . || { echo 's3cmd get failed' ; exit 1; }
if [ -n "${3}" ] ; then
    mv $fn ${3} || true
fi
"""
                )
        base.copy_bootstrap = path_join(
                                        True, base.dependency_dir,
                                        os.path.basename(
                                                copy_bootstrap
                                            )
                                    )
        ansible.put(copy_bootstrap, base.copy_bootstrap)
        os.remove(copy_bootstrap)
        if base.secure_stack_name is not None:
            encrypt_bootstrap = os.path.join(temp_dependency_dir,
                                                'encrypt.sh')
            with open(encrypt_bootstrap, 'w') as script_stream:
                # Code taken from http://bddubois-emr.s3.amazonaws.com/emr-volume-encryption.sh
                print >>script_stream, (
"""#!/usr/bin/env bash
set -ex

EPHEMERAL_MNT_DIRS=`awk '/mnt/{print $2}' < /proc/mounts`
ENCRYPTED_SIZE=8g
export PASSWORD=$(dd count=3 bs=16 if=/dev/urandom of=/dev/stdout 2>/dev/null | base64) PASSWORD_FILE="/tmp/pwd"
i=0
STATUS=0
TMPSIZE=40

echo ${PASSWORD} > $PASSWORD_FILE
if [ ! $? -eq 0 ]; then
  echo "ERROR: Failed to create password file"
  STATUS=1
fi

sudo modprobe loop

# Install cryptsetup
#THESE steps fail but its ok the clone in amazon ami already has this installed # #sudo yum -y update #sudo yum -y install cryptsetup

mychildren=""

if [ $STATUS -eq 0 ]; then
  for DIR in $EPHEMERAL_MNT_DIRS; do
    #
    # Set up some variables
    #
    ENCRYPTED_LOOPBACK_DIR=$DIR/encrypted_loopbacks
    ENCRYPTED_MOUNT_POINT=$DIR/space.encrypted/
    ENCRYPTED_SPACE=$DIR/space
    DFS_DATA_DIR=$DIR/var/lib/hadoop/dfs
    TMP_DATA_DIR=$DIR/var/lib/hadoop/tmp
    ENCRYPTED_LOOPBACK_DEVICE=/dev/loop$i
    ENCRYPTED_NAME=crypt$i

    mkdir -p ${ENCRYPTED_SPACE}
    
    if [ $STATUS -eq 0 ]; then
      # Get the total number of blocks for this filesystem $DIR
      nblocks=`stat -f -c '%a' $DIR`
      # Get the block size (in bytes for this filesystem $DIR)
      bsize=`stat -f -c '%s' $DIR`
      # Calculate the mntsize in MB (divisible by 1000)
      mntsize=`expr $nblocks \* $bsize \/ 1000 \/ 1000 \/ 1000`
      # Make $TMPSIZE 1/10 of mntsize
      TMPSIZE=`expr $mntsize \/ 10`
      if [ ! $? -eq 0 ]; then
        echo "ERROR: Failed to get mount size"
        STATUS=1
      fi
      ENCRYPTED_SIZE=`expr $mntsize - $TMPSIZE`g
      if [ ! $? -eq 0 ]; then
        echo "ERROR: Failed to calculate encrypted size"
        STATUS=1
      fi
    fi

    #
    # Create directories
    #
    if [ $STATUS -eq 0 ]; then
      mkdir $ENCRYPTED_LOOPBACK_DIR
      if [ ! $? -eq 0 ]; then
        echo "ERROR: Failed to get create ENCRYPTED_LOOPBACK_DIR"
        STATUS=1
      else 
        echo SUCCESS: Created directory $ENCRYPTED_LOOPBACK_DIR
      fi
    fi
    if [ $STATUS -eq 0 ]; then
      mkdir $ENCRYPTED_MOUNT_POINT
      if [ ! $? -eq 0 ]; then
        echo "ERROR: Failed to get create ENCRYPTED_MOUNT_POINT"
        STATUS=1
      else
        echo SUCCESS: Created directory $ENCRYPTED_MOUNT_POINT
      fi
    fi
    #
    # Create loopback device
    #
    if [ $STATUS -eq 0 ]; then
      sudo fallocate -l $ENCRYPTED_SIZE $ENCRYPTED_LOOPBACK_DIR/encrypted_loopback.img
      if [ ! $? -eq 0 ]; then
        echo "ERROR: Failed to allocate $ENCRYPTED_SIZE $ENCRYPTED_LOOPBACK_DIR/encrypted_loopback.img"
        STATUS=1
      else
        echo SUCCESS: Allocated $ENCRYPTED_SIZE $ENCRYPTED_LOOPBACK_DIR/encrypted_loopback.img
      fi
    fi 
    if [ $STATUS -eq 0 ]; then
      sudo chown hadoop:hadoop $ENCRYPTED_LOOPBACK_DIR/encrypted_loopback.img
      if [ ! $? -eq 0 ]; then
        echo "ERROR: Failed to chown $ENCRYPTED_LOOPBACK_DIR/encrypted_loopback.img"
        STATUS=1
      else
        echo SUCCESS: chowned $ENCRYPTED_LOOPBACK_DIR/encrypted_loopback.img
      fi
    fi
    if [ $STATUS -eq 0 ]; then
      sudo losetup /dev/loop$i $ENCRYPTED_LOOPBACK_DIR/encrypted_loopback.img
      if [ ! $? -eq 0 ]; then
        echo "ERROR: Failed to losetup $ENCRYPTED_LOOPBACK_DIR/encrypted_loopback.img"
        STATUS=1
      else
        echo SUCCESS: losetup $ENCRYPTED_LOOPBACK_DIR/encrypted_loopback.img
      fi
    fi
    #
    # Set up LUKS
    #
    if [ $STATUS -eq 0 ]; then
      sudo cryptsetup luksFormat -q --key-file $PASSWORD_FILE $ENCRYPTED_LOOPBACK_DEVICE
      if [ ! $? -eq 0 ]; then
        echo "ERROR: Failed to  cryptsetup luksFormat -q --key-file $PASSWORD_FILE $ENCRYPTED_LOOPBACK_DEVICE"
        STATUS=1
      else
        echo SUCCESS: cryptsetup luksFormat 
      fi
    fi
    if [ $STATUS -eq 0 ]; then
      sudo cryptsetup luksOpen -q --key-file $PASSWORD_FILE $ENCRYPTED_LOOPBACK_DEVICE $ENCRYPTED_NAME
      if [ ! $? -eq 0 ]; then
        echo "ERROR: Failed to  cryptsetup luksOpen -q --key-file $PASSWORD_FILE $ENCRYPTED_LOOPBACK_DEVICE"
        STATUS=1
      else
        echo SUCCESS: cryptsetup luksopen
      fi
    fi

    #
    # Create file system
    #
    if [ $STATUS -eq 0 ]; then
    mycmd="sudo mkfs.ext4 -m 0 -E lazy_itable_init=1 /dev/mapper/$ENCRYPTED_NAME && sudo mount /dev/mapper/$ENCRYPTED_NAME $ENCRYPTED_SPACE && sudo mkdir -p $ENCRYPTED_SPACE/dfs && sudo mkdir -p $ENCRYPTED_SPACE/tmp/nm-local-dir && sudo rm -rf $DFS_DATA_DIR && sudo ln -s $ENCRYPTED_SPACE/dfs $DFS_DATA_DIR && sudo chown hadoop:hadoop $ENCRYPTED_SPACE/dfs && sudo chown hadoop:hadoop $DFS_DATA_DIR && sudo rm -rf $DFS_DATA_DIR/lost\+found && sudo rm -rf $TMP_DATA_DIR && sudo ln -s $ENCRYPTED_SPACE/tmp $TMP_DATA_DIR && sudo chown hadoop:hadoop $ENCRYPTED_SPACE/tmp && sudo chown hadoop:hadoop $TMP_DATA_DIR && sudo chown hadoop:hadoop $TMP_DATA_DIR/nm-local-dir && sudo echo iamdone-$ENCRYPTED_NAME && sudo chown -R hadoop:hadoop $ENCRYPTED_SPACE && sudo chmod -R 0755 $ENCRYPTED_SPACE && date "
    echo $mycmd
    eval $mycmd &
      if [ ! $? -eq 0 ]; then
        echo "ERROR: Failed to run the my cmd that follows $mycmd"
        STATUS=1
      else
        echo SUCCESS: MYCMD 
      fi
    fi

    mychildren="$mychildren $!"

    let i=i+1
done
fi

for mypid in $mychildren
do
    wait $mypid
done

sudo rm -f $PASSWORD_FILE

date
echo "everything done"
echo $STATUS STATUS
exit $STATUS
"""
            )
            base.encrypt_bootstrap = path_join(
                                            True, base.dependency_dir,
                                            os.path.basename(
                                                    encrypt_bootstrap
                                                )
                                        )
            ansible.put(encrypt_bootstrap, base.encrypt_bootstrap)
            os.remove(encrypt_bootstrap)
        if sra_tools_needed:
            vdb_bootstrap = os.path.join(temp_dependency_dir,
                                                'vdb.sh')
            if base.dbgap_present:
                with open(vdb_bootstrap, 'w') as script_stream:
                    print >>script_stream, (
"""#!/usr/bin/env bash

mkdir -p {vdb_workspace}/insecure
mkdir -p ~/.ncbi
cat >~/.ncbi/user-settings.mkfg <<EOF
/repository/user/main/public/root = "{vdb_workspace}/insecure"
EOF
mkdir -p {vdb_workspace}/secure
{vdb_config} --import /mnt/space/DBGAP.ngc {vdb_workspace}/secure
"""
                    ).format(vdb_config=_elastic_vdb_config_exe,
                             vdb_workspace=_elastic_vdb_workspace)
            else:
                # Don't need secure dir for workspace
                with open(vdb_bootstrap, 'w') as script_stream:
                    print >>script_stream, (
"""#!/usr/bin/env bash

mkdir -p {vdb_workspace}/insecure
mkdir -p ~/.ncbi
cat >~/.ncbi/user-settings.mkfg <<EOF
/repository/user/main/public/root = "{vdb_workspace}/insecure"
EOF
"""
                    ).format(vdb_workspace=_elastic_vdb_workspace)
        base.vdb_bootstrap = path_join(
                                        True, base.dependency_dir,
                                        os.path.basename(
                                                vdb_bootstrap
                                            )
                                    )
        ansible.put(vdb_bootstrap, base.vdb_bootstrap)
        os.remove(vdb_bootstrap)
        shutil.rmtree(temp_dependency_dir)
        print_to_screen('Copied Rail-RNA and bootstraps to S3.',
                         newline=True, carriage_return=False)
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
        if not (float(master_instance_count).is_integer()
                and master_instance_count >= 1):
            base.errors.append('Master instance count '
                               '(--master-instance-count) must be an '
                               'integer >= 1, but {0} was entered.'.format(
                                                    master_instance_count
                                                ))
        base.master_instance_count = master_instance_count
        if not (float(core_instance_count).is_integer()
                 and core_instance_count >= 1):
            base.errors.append('Core instance count '
                               '(--core-instance-count) must be an '
                               'integer >= 1, but {0} was entered.'.format(
                                                    core_instance_count
                                                ))
        base.core_instance_count = core_instance_count
        if not (float(task_instance_count).is_integer()
                and task_instance_count >= 0):
            base.errors.append('Task instance count '
                               '(--task-instance-count) must be an '
                               'integer >= 1, but {0} was entered.'.format(
                                                    task_instance_count
                                                ))
        base.task_instance_count = task_instance_count
        # Raise exceptions before computing mems
        raise_runtime_error(base)
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
        base.consistent_view = consistent_view
        base.no_direct_copy = no_direct_copy

    @staticmethod
    def add_args(general_parser, required_parser, output_parser, 
                    elastic_parser, align=False):
        if align:
            required_parser.add_argument(
                '-i', '--input', type=str, required=True, metavar='<s3_dir>',
                help='input directory with preprocessed reads; must begin ' \
                     'with s3://'
            )
        else:
            # "prep" or "go" flows
            general_parser.add_argument(
                '--dbgap-key', required=False, metavar='<file>',
                default=None,
                help='path to dbGaP key file, which has the extension "ngc" '
                     'for SRA data and "key" for CGHub data; '
                     'must be on local filesystem. This file is uploaded '
                     'securely to S3 and scheduled for deletion if dbGaP '
                     'data is present in the manifest file'
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
        general_parser.add_argument(
            '--secure-stack-name', type=str, required=False,
            metavar='<str>',
            default=None,
            help='name of secure stack enforcing dbGaP requirements; '
                 'must be specified when analyzing dbGaP data.'
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
            default='3.8.0',
            help='Amazon Machine Image version'
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
        elastic_parser.add_argument('--consistent-view',
            action='store_const',
            const=True,
            default=False,
            help=('use "consistent view," which incurs DynamoDB '
                  'charges; not necessary for Rail\'s default job flows'))
        elastic_parser.add_argument('--no-direct-copy',
            action='store_const',
            const=True,
            default=False,
            help=('writes intermediate data to HDFS before copying to S3; '
                  'helps ensure that low-probability data loss between steps '
                  'does not occur'))
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
                  'should be spot (def: none)')
        )
        elastic_parser.add_argument('--core-instance-bid-price', type=float,
            metavar='<dec>',
            required=False,
            default=None,
            help=('bid price (dollars/hr); invoke only if core instances '
                  'should be spot (def: none)')
        )
        elastic_parser.add_argument('--task-instance-bid-price', type=float,
            metavar='<dec>',
            required=False,
            default=None,
            help=('bid price (dollars/hr); invoke only if task instances '
                  'should be spot (def: none)')
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
            help=('core instance type (def: master instance type)')
        )
        elastic_parser.add_argument('--task-instance-type', type=str,
            metavar='<choice>',
            required=False,
            default=None,
            help=('task instance type (def: master instance type)')
        )
        elastic_parser.add_argument('--ec2-key-name', type=str,
            metavar='<str>',
            required=False,
            default=None,
            help=('key pair name for SSHing to EC2 instances (def: '
                  'unspecified, so SSHing is not permitted)')
        )
        elastic_parser.add_argument('--ec2-subnet-id', type=str,
            metavar='<str>',
            required=False,
            default=None,
            help=('ID of subnet into which EMR cluster should be '
                  'launched; a properly configured VPC subnet '
                  'must be specified in secure mode (def: none)'))
        elastic_parser.add_argument('--ec2-master-security-group-id', type=str,
            metavar='<str>',
            required=False,
            default=None,
            help=('ID of security group to associate with every master '
                  'instance in the EMR cluster (def: Amazon default)'))
        elastic_parser.add_argument('--ec2-slave-security-group-id', type=str,
            metavar='<str>',
            required=False,
            default=None,
            help=('ID of security group to associate with every slave '
                  'instance in the EMR cluster (def: Amazon default)'))
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
            default=None,
            help=('Amazon data center in which to run job flow. Google '
                  '"Elastic MapReduce regions" for recent list of centers '
                  '(def: from --profile)')
        )
        elastic_parser.add_argument('--service-role', type=str,
            metavar='<str>',
            required=False,
            default=None,
            help=('IAM service role (def: from --profile if available; '
                  'otherwise, attempts "EMR_DefaultRole")')
        )
        elastic_parser.add_argument('--instance-profile', type=str,
            metavar='<str>',
            required=False,
            default=None,
            help=('IAM EC2 instance profile (def: from --profile if available'
                  '; otherwise, attempts "EMR_EC2_DefaultRole")')
        )
        general_parser.add_argument(
            '--max-task-attempts', type=int, required=False,
            metavar='<int>',
            default=4,
            help=('maximum number of attempts per task; sets both '
                  'mapreduce.map.maxattempts and mapreduce.reduce.maxattempts')
        )

    @staticmethod
    def hadoop_debugging_steps(base):
        return [
            {
                'ActionOnFailure' : base.action_on_failure,
                'HadoopJarStep' : {
                    'Args' : [
                        ('s3://%s.elasticmapreduce/libs/'
                         'state-pusher/0.1/fetch') % base.region
                    ],
                    'Jar' : ('s3://%s.elasticmapreduce/libs/'
                             'script-runner/script-runner.jar') % base.region
                },
                'Name' : 'Set up Hadoop Debugging'
            }
        ]

    @staticmethod
    def prebootstrap(base):
        if base.secure_stack_name is not None:
            to_return = [
                    {
                        'Name' : 'Encrypt ephemeral space',
                        'ScriptBootstrapAction' : {
                            'Args' : [],
                            'Path' : base.encrypt_bootstrap
                        }
                    }
                ]
        else:
            to_return = []
        return to_return + [
                {
                    'Name' : 'Install S3cmd',
                    'ScriptBootstrapAction' : {
                        'Args' : [],
                        'Path' : base.install_s3cmd_bootstrap
                    }
                }
            ]

    @staticmethod
    def bootstrap(base):
        return [
            {
                'Name' : 'Allocate swap space',
                'ScriptBootstrapAction' : {
                    'Args' : [
                        '%d' % base.mem,
                        '/mnt/space/swapfile'
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
                        % (base.nodemanager_mem / base.max_tasks),
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
                        % (base.nodemanager_mem / base.max_tasks),
                        '-m',
                        'mapreduce.reduce.memory.mb=%d'
                        % (base.nodemanager_mem / base.max_tasks),
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
                        'mapreduce.reduce.maxattempts=%d' 
                        % base.max_task_attempts,
                        '-m',
                        'mapreduce.map.maxattempts=%d'
                        % base.max_task_attempts,
                        '-m',
                        ('mapreduce.output.fileoutputformat.compress.codec='
                         'com.hadoop.compression.lzo.LzopCodec'),
                        '-m',
                        'mapreduce.job.maps=%d' % base.total_cores,
                        '-e',
                        'fs.s3.enableServerSideEncryption=true'
                    ] + (['-e', 'fs.s3.consistent=true']
                            if base.consistent_view
                            else ['-e', 'fs.s3.consistent=false']),
                    'Path' : ('s3://%s.elasticmapreduce/bootstrap-actions/'
                              'configure-hadoop' % base.region)
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
        if base.ec2_subnet_id is not None:
            to_return['Ec2SubnetId'] = base.ec2_subnet_id
        if base.ec2_master_security_group_id is not None:
            to_return['EmrManagedMasterSecurityGroup'] \
                = base.ec2_master_security_group_id
        if base.ec2_slave_security_group_id is not None:
            to_return['EmrManagedSlaveSecurityGroup'] \
                = base.ec2_slave_security_group_id
        return to_return

class RailRnaPreprocess(object):
    """ Sets parameters relevant to just the preprocessing step of a job flow.
    """
    def __init__(self, base, nucleotides_per_input=8000000, gzip_input=True,
                    do_not_bin_quals=False, short_read_names=False,
                    skip_bad_records=False):
        if not (float(nucleotides_per_input).is_integer() and
                nucleotides_per_input > 0):
            base.errors.append('Nucleotides per input '
                               '(--nucleotides-per-input) must be an integer '
                               '> 0, but {0} was entered.'.format(
                                                        nucleotides_per_input
                                                       ))
        base.nucleotides_per_input = nucleotides_per_input
        base.gzip_input = gzip_input
        base.do_not_bin_quals = do_not_bin_quals
        base.skip_bad_records = skip_bad_records
        base.short_read_names = short_read_names


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
        output_parser.add_argument(
                '--do-not-bin-quals', action='store_const', const=True,
                default=False,
                help=('does not place phred quality score of each read base '
                      'into one of five bins; binning saves space without '
                      'affecting alignment quality')
            )
        output_parser.add_argument(
                '--short-read-names', action='store_const', const=True,
                default=False,
                help=('use short read names to save space; this loses '
                      'original read names but preserves information about '
                      'which reads are pairs')
            )
        output_parser.add_argument(
                '--skip-bad-records', action='store_const', const=True,
                default=False,
                help=('skips all bad records rather than raising an exception '
                      'on encountering a bad record')
            )
        general_parser.add_argument(
            '--do-not-check-manifest', action='store_const', const=True,
            default=False,
            help='do not check that files listed in manifest file exist'
        )

    @staticmethod
    def protosteps(base, prep_dir, push_dir, elastic=False):
        if not elastic:
            steps_to_return = [
                {
                    'name' : 'Count lines in input files',
                    'mapper' : 'count_inputs.py',
                    'inputs' : [base.old_manifest
                                if hasattr(base, 'old_manifest')
                                else base.manifest],
                    'no_input_prefix' : True,
                    'output' : 'count_lines',
                    'inputformat' : (
                           'org.apache.hadoop.mapred.lib.NLineInputFormat'
                        )
                },
                {
                    'name' : 'Assign reads to preprocessing tasks',
                    'reducer' : ('assign_splits.py --num-processes {0}'
                                 ' --out {1} --filename {2} {3}').format(
                                                        base.num_processes,
                                                        base.intermediate_dir,
                                                        'split.manifest',
                                                        ('--scratch %s' %
                                                          base.scratch)
                                                        if base.scratch
                                                        is not None
                                                        else ''
                                                    ),
                    'inputs' : ['count_lines'],
                    'output' : 'assign_reads',
                    'tasks' : 1,
                    'partition' : '-k1,1'
                },
                {
                    'name' : 'Preprocess reads',
                    'mapper' : ('preprocess.py --nucs-per-file={0} {1} '
                                '--push={2} --gzip-level {3} {4} {5} '
                                '{6} {7} {8} {9}').format(
                                                base.nucleotides_per_input,
                                                '--gzip-output' if
                                                base.gzip_input else '',
                                                push_dir,
                                                base.gzip_level if
                                                'gzip_level' in
                                                dir(base) else 3,
                                                '--stdout' if elastic
                                                else '',
                                                ('--scratch %s' %
                                                  base.scratch) if 
                                                base.scratch is not None
                                                else '',
                                                '--bin-qualities'
                                                if
                                                not base.do_not_bin_quals
                                                else '',
                                                '--shorten-read-names'
                                                if base.short_read_names
                                                else '',
                                                '--skip-bad-records'
                                                if base.skip_bad_records
                                                else '',
                                                '--fastq-dump-exe %s'
                                                    % (base.fastq_dump_exe
                                                        if hasattr(base,
                                                            'fastq_dump_exe')
                                                        else 'fastq-dump')
                                            ),
                    'inputs' : [os.path.join(base.intermediate_dir,
                                                'split.manifest')],
                    'no_input_prefix' : True,
                    'output' : push_dir if elastic else prep_dir,
                    'no_output_prefix' : True,
                    'inputformat' : (
                           'org.apache.hadoop.mapred.lib.NLineInputFormat'
                        )
                },
            ]
        else:
            steps_to_return = [
                {
                    'name' : 'Preprocess reads',
                    'mapper' : ('preprocess.py --nucs-per-file={0} {1} '
                                '--push={2} --gzip-level {3} {4} {5} {6} {7} '
                                '--keep-alive {8} {9}').format(
                                                base.nucleotides_per_input,
                                                '--gzip-output' if
                                                base.gzip_input else '',
                                                ab.Url(push_dir).to_url(
                                                        caps=True
                                                    ),
                                                base.gzip_level if
                                                'gzip_level' in
                                                dir(base) else 3,
                                                '--stdout' if elastic
                                                else '',
                                                '--bin-qualities'
                                                if
                                                not base.do_not_bin_quals
                                                else '',
                                                '--shorten-read-names'
                                                if base.short_read_names
                                                else '',
                                                '--skip-bad-records'
                                                if base.skip_bad_records
                                                else '',
                                            ('--workspace-dir %s/secure'
                                             % _elastic_vdb_workspace)
                                                if base.dbgap_present
                                                else '',
                                                '--fastq-dump-exe %s'
                                                    % (base.fastq_dump_exe
                                                        if hasattr(base,
                                                            'fastq_dump_exe')
                                                        else 'fastq-dump')
                                            ),
                    'inputs' : [base.old_manifest
                                if hasattr(base, 'old_manifest')
                                else base.manifest],
                    'no_input_prefix' : True,
                    'output' : push_dir if elastic else prep_dir,
                    'no_output_prefix' : True,
                    'inputformat' : (
                           'org.apache.hadoop.mapred.lib.NLineInputFormat'
                        )
                },
            ]
        return steps_to_return

    @staticmethod
    def srabootstrap(base):
        return [
                    {
                        'Name' : 'Set up SRA Tools workspace',
                        'ScriptBootstrapAction' : {
                            'Args' : [],
                            'Path' : base.vdb_bootstrap
                        }
                }
            ]

    @staticmethod
    def bootstrap(base):
        return [
            {
                'Name' : 'Install Rail-RNA and create JAR dependencies',
                'ScriptBootstrapAction' : {
                    'Args' : [
                        base.elastic_rail_path, _jar_target,
                        '-y', '-p', '-s'
                    ], # always say yes, only prep dependencies, and symlink  
                    'Path' : base.install_rail_bootstrap
                }
            }
        ]

class RailRnaAlign(object):
    """ Sets parameters relevant to just the "align" job flow. """
    def __init__(self, base, input_dir=None, elastic=False,
        bowtie1_exe=None, bowtie_idx='genome', bowtie1_build_exe=None,
        bowtie2_exe=None, bowtie2_build_exe=None, k=1, bowtie2_args='',
        samtools_exe=None, bedgraphtobigwig_exe=None,
        partition_length=5000, max_readlet_size=25,
        readlet_config_size=32, min_readlet_size=15, readlet_interval=4,
        cap_size_multiplier=1.2, max_intron_size=500000, min_intron_size=10,
        min_exon_size=9, library_size=40, search_filter='none',
        motif_search_window_size=1000, max_gaps_mismatches=3,
        motif_radius=5, genome_bowtie1_args='-v 0 -a -m 80',
        transcriptome_bowtie2_args='-k 30', count_multiplier=15,
        junction_criteria='0.5,5', indel_criteria='0.5,5', tie_margin=6,
        normalize_percentile=0.75, transcriptome_indexes_per_sample=500,
        drop_deletions=False, do_not_output_bam_by_chr=False,
        do_not_output_ave_bw_by_chr=False, output_sam=False,
        do_not_drop_polyA_tails=False, deliverables='idx,tsv,bed,bw',
        bam_basename='alignments', bed_basename='', tsv_basename='',
        assembly='hg19', s3_ansible=None):
        if not elastic:
            '''Programs and Bowtie indexes should be checked only in local
            mode. First grab Bowtie index paths.'''
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
            bowtie_idx = ','.join(bowtie_idx)
            if ',' in bowtie_idx:
                bowtie1_idx, _, bowtie2_idx = bowtie_idx.partition(',')
            else:
                bowtie1_idx = bowtie2_idx = bowtie_idx
            bowtie1_idx = os.path.expandvars(os.path.expanduser(bowtie1_idx))
            bowtie2_idx = os.path.expandvars(os.path.expanduser(bowtie2_idx))
            bowtie1_extensions = set(['.1.ebwt', '.2.ebwt', '.3.ebwt',
                                      '.4.ebwt', '.rev.1.ebwt', '.rev.2.ebwt'])
            bowtie2_extensions = set(['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', 
                                      '.rev.1.bt2', '.rev.2.bt2'])
            if not (ab.Url(bowtie1_idx).is_local
                    and ab.Url(bowtie2_idx).is_local):
                base.errors.append(('Bowtie index basenames must be '
                                    'on the local filesystem, but "{0}" and '
                                    '"{1}" were entered.').format(
                                            bowtie1_idx, bowtie2_idx
                                        ))
            else:
                import glob
                bowtie1_idx_len = len(bowtie1_idx)
                existing_bowtie1_extensions = set(
                        [idx_file[bowtie1_idx_len:]
                            for idx_file in glob.glob(bowtie1_idx + '*')]
                    )
                if not (bowtie1_extensions <= existing_bowtie1_extensions):
                    # Swap BT1 and 2 indexes and try again
                    bowtie1_idx, bowtie2_idx = bowtie2_idx, bowtie1_idx
                bowtie1_idx_len = len(bowtie1_idx)
                existing_bowtie1_extensions = set(
                        [idx_file[bowtie1_idx_len:]
                            for idx_file in glob.glob(bowtie1_idx + '*ebwt')]
                    )
                bowtie2_idx_len = len(bowtie2_idx)
                existing_bowtie2_extensions = set(
                        [idx_file[bowtie2_idx_len:]
                            for idx_file in glob.glob(bowtie2_idx + '*bt2')]
                    )
                missing_extensions = ['"{}"'.format(extension)
                                            for extension in
                                            ((bowtie1_extensions
                                                | bowtie2_extensions) - (
                                                existing_bowtie2_extensions
                                                | existing_bowtie1_extensions
                                            ))]
                if missing_extensions:
                    base.errors.append(('Some Bowtie index files were not '
                                        'found. Check that index files with '
                                        'the extensions {0} exist with the '
                                        'basename(s) specified.').format(
                                                ', '.join(
                                                    missing_extensions[:-1]
                                                ) + ', and '
                                                + missing_extensions[-1]
                                            ))
            base.bowtie1_idx, base.bowtie2_idx = bowtie1_idx, bowtie2_idx
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
            # Output any errors before detect message is determined
            raise_runtime_error(base)
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
            try:
                base.index_archive = _assemblies[assembly]
            except KeyError:
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
            # Don't check these dependencies
            base.bedgraphtobigwig_exe = _elastic_bedgraphtobigwig_exe
            base.samtools_exe = _elastic_samtools_exe
            base.bowtie1_exe = _elastic_bowtie1_exe
            base.bowtie2_exe = _elastic_bowtie2_exe
            base.bowtie1_build_exe = _elastic_bowtie1_build_exe
            base.bowtie2_build_exe = _elastic_bowtie2_build_exe
        if '--mp' in bowtie2_args:
            print_to_screen('Warning: --mp parameter in specified in Bowtie 2 '
                            'arguments (--bowtie2-args) will be ignored.',
                            newline=True, carriage_return=False)
            # Replace -k value in bowtie2 args with -k arg
            bowtie2_arg_tokens = [
                    token.strip() for token
                    in bowtie2_args.replace('=', ' ').split(' ')
                ]
            # Remove --mp parameter
            bowtie2_args = ' '.join(
                        [token for i, token in enumerate(bowtie2_arg_tokens)
                         if (token != '--mp'
                             and bowtie2_arg_tokens[max(i-1, 0)] != '--mp')]
                    )
        if k is not None:
            if not (float(k).is_integer() and k >= 1):
                base.errors.append('Number of alignments to report per read '
                                   '(-k) must be an integer >= 1, but '
                                   '{0} was entered.'.format(k))
            else:
                # Replace -k value in bowtie2 args with -k arg
                bowtie2_arg_tokens = [
                        token.strip() for token
                        in bowtie2_args.replace('=', ' ').split(' ')
                    ]
                try:
                    bowtie2_arg_tokens[bowtie2_arg_tokens.index('-k') + 1] \
                        = str(k)
                except ValueError:
                    if bowtie2_arg_tokens[0] == '':
                        bowtie2_arg_tokens = ['-k', str(k)]
                    else:
                        bowtie2_arg_tokens.extend(['-k', str(k)])
                except IndexError:
                    bowtie2_arg_tokens.append(str(k))
                base.bowtie2_args = ' '.join(bowtie2_arg_tokens)
        else:
            base.bowtie2_args = bowtie2_args
        base.k = k
        if not (float(partition_length).is_integer() and
                partition_length > 0):
            base.errors.append('Genome partition length '
                               '(--partition-length) must be an '
                               'integer > 0, but {0} was entered.'.format(
                                                        partition_length
                                                    ))
        base.partition_length = partition_length
        if not (float(min_readlet_size).is_integer() and min_readlet_size > 0):
            base.errors.append('Minimum readlet size (--min-readlet-size) '
                               'must be an integer > 0, but '
                               '{0} was entered.'.format(min_readlet_size))
        base.min_readlet_size = min_readlet_size
        if not (float(max_readlet_size).is_integer() and max_readlet_size
                >= min_readlet_size):
            base.errors.append('Maximum readlet size (--max-readlet-size) '
                               'must be an integer >= minimum readlet size '
                               '(--min-readlet-size) = '
                               '{0}, but {1} was entered.'.format(
                                                    base.min_readlet_size,
                                                    max_readlet_size
                                                ))
        base.max_readlet_size = max_readlet_size
        if not (float(readlet_config_size).is_integer() and readlet_config_size
                >= max_readlet_size):
            base.errors.append('Readlet config size (--readlet-config-size) '
                               'must be an integer >= maximum readlet size '
                               '(--max-readlet-size) = '
                               '{0}, but {1} was entered.'.format(
                                                    base.max_readlet_size,
                                                    readlet_config_size
                                                ))
        base.readlet_config_size = readlet_config_size
        if not (float(readlet_interval).is_integer() and readlet_interval > 0):
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
        if not (float(min_intron_size).is_integer() and min_intron_size > 0):
            base.errors.append('Minimum intron size (--min-intron-size) '
                               'must be an integer > 0, but '
                               '{0} was entered.'.format(min_intron_size))
        base.min_intron_size = min_intron_size
        if not (float(max_intron_size).is_integer() and max_intron_size
                >= min_intron_size):
            base.errors.append('Maximum intron size (--max-intron-size) '
                               'must be an integer >= minimum intron size '
                               '(--min-readlet-size) = '
                               '{0}, but {1} was entered.'.format(
                                                    base.min_intron_size,
                                                    max_intron_size
                                                ))
        base.max_intron_size = max_intron_size
        if not (float(min_exon_size).is_integer() and min_exon_size > 0):
            base.errors.append('Minimum exon size (--min-exon-size) '
                               'must be an integer > 0, but '
                               '{0} was entered.'.format(min_exon_size))
        base.min_exon_size = min_exon_size
        if search_filter == 'none':
            base.search_filter = 1
        elif search_filter == 'mild':
            base.search_filter = int(base.min_exon_size * 2 / 3)
        elif search_filter == 'strict':
            try:
                base.search_filter = base.min_exon_size
            except:
                pass
        else:
            try:
                base.search_filter = int(search_filter)
            except ValueError:
                # Not an integer
                base.errors.append('Search filter (--search-filter) '
                                   'must be an integer >= 1 or one of '
                                   '{"none", "mild", "strict"}, but {0} was '
                                   'entered.'.format(search_filter))
        if not (float(motif_search_window_size).is_integer() and 
                    motif_search_window_size >= 0):
            base.errors.append('Motif search window size '
                               '(--motif-search-window-size) must be an '
                               'integer >= 0, but {0} was entered.'.format(
                                                    motif_search_window_size
                                                ))
        base.motif_search_window_size = motif_search_window_size
        if max_gaps_mismatches is not None and not (
                float(max_gaps_mismatches).is_integer() and 
                max_gaps_mismatches >= 0
            ):
            base.errors.append('Max gaps and mismatches '
                               '(--max-gaps-mismatches) must be an '
                               'integer >= 0, but {0} was entered.'.format(
                                                    max_gaps_mismatches
                                                ))
        base.max_gaps_mismatches = max_gaps_mismatches
        if not (float(motif_radius).is_integer() and
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
        if not (float(tie_margin).is_integer() and
                    tie_margin >= 0):
            base.errors.append('Tie margin (--tie-margin) must be an '
                               'integer >= 0, but {0} was entered.'.format(
                                                    tie_margin
                                                ))
        base.tie_margin = tie_margin
        if not (float(count_multiplier).is_integer() and
                    count_multiplier >= 0):
            base.errors.append('Count multiplier (--count-multiplier) must '
                               'be an integer >= 0, but '
                               '{0} was entered.'.format(
                                                    count_multiplier
                                                ))
        base.count_multiplier = count_multiplier
        if not (float(library_size).is_integer() and
                    library_size >= 0):
            base.errors.append('Library size in millions of reads '
                               '(--library-size) must be an integer >= 0, but '
                               '{0} was entered.'.format(
                                                    library_size
                                                ))
        base.library_size = library_size
        if isinstance(junction_criteria, str):
            junction_criteria = [junction_criteria]
        confidence_criteria_split = [criterion.strip(_whitespace_and_comma)
                                        for criterion in 
                                        ','.join(junction_criteria).split(',')
                                        if criterion]
        confidence_criteria_error = len(confidence_criteria_split) != 2
        if not confidence_criteria_error:
            try:
                base.sample_fraction = float(confidence_criteria_split[0])
            except ValueError:
                confidence_criteria_error = True
            else:
                if not (0 <= base.sample_fraction <= 1):
                    confidence_criteria_error = True
            try:
                base.coverage_threshold = int(confidence_criteria_split[1])
            except ValueError:
                confidence_criteria_error = True
            else:
                if not (base.coverage_threshold >= 0
                        or base.coverage_threshold == -1):
                    confidence_criteria_error = True
        if confidence_criteria_error:
            base.errors.append('Junction confidence criteria '
                               '(--junction-criteria) must be a '
                               'comma- or space-separated list of two '
                               'elements: the first should be a decimal value '
                               'between 0 and 1 inclusive; the second should '
                               'be an integer >= 1 or -1 to disable filtering '
                               'by read count. {0} was entered.'.format(
                                                    ','.join(junction_criteria)
                                                ))
        if isinstance(indel_criteria, str):
            indel_criteria = [indel_criteria]
        confidence_criteria_split = [criterion.strip(_whitespace_and_comma)
                                        for criterion in 
                                        ','.join(indel_criteria).split(',')
                                        if criterion]
        confidence_criteria_error = len(confidence_criteria_split) != 2
        if not confidence_criteria_error:
            try:
                base.indel_sample_fraction = float(
                                            confidence_criteria_split[0]
                                        )
            except ValueError:
                confidence_criteria_error = True
            else:
                if not (0 <= base.indel_sample_fraction <= 1):
                    confidence_criteria_error = True
            try:
                base.indel_coverage_threshold = int(
                                            confidence_criteria_split[1]
                                        )
            except ValueError:
                confidence_criteria_error = True
            else:
                if not (base.indel_coverage_threshold >= 0
                        or base.indel_coverage_threshold == -1):
                    confidence_criteria_error = True
        if confidence_criteria_error:
            base.errors.append('Indel confidence criteria '
                               '(--indel-criteria) must be a '
                               'comma- or space-separated list of two '
                               'elements: the first should be a decimal value '
                               'between 0 and 1 inclusive; the second should '
                               'be an integer >= 1 or -1 to disable filtering '
                               'by read count. {0} was entered.'.format(
                                                    ','.join(indel_criteria)
                                                ))
        base.do_not_drop_polyA_tails = do_not_drop_polyA_tails
        base.drop_deletions = drop_deletions
        base.do_not_output_bam_by_chr = do_not_output_bam_by_chr
        base.do_not_output_ave_bw_by_chr = do_not_output_ave_bw_by_chr
        if not (float(transcriptome_indexes_per_sample).is_integer()
                    and (1000 >= transcriptome_indexes_per_sample >= 1)):
            base.errors.append('Transcriptome indexes per sample '
                               '(--transcriptome-indexes-per-sample) must be '
                               'an integer between 1 and 1000.')
        base.transcriptome_indexes_per_sample \
            = transcriptome_indexes_per_sample
        base.bam_basename = bam_basename
        base.bed_basename = bed_basename
        base.tsv_basename = tsv_basename
        deliverable_choices = set(
                ['idx', 'bam', 'sam', 'bed', 'tsv', 'bw', 'jx']
            )
        if isinstance(deliverables, str):
            deliverables = [deliverables]
        split_deliverables = set([deliverable.strip(_whitespace_and_comma)
                                    for deliverable
                                    in ','.join(deliverables).split(',')
                                    if deliverable])
        undeliverables = split_deliverables - deliverable_choices
        if undeliverables:
            base.errors.append('Some deliverables (--deliverables) specified '
                               'are invalid. Valid choices are in {{"idx", '
                               '"bam", "bed", "tsv", "bw", "jx"}}, but '
                               '"{0}" was entered.'.format(deliverables))
        elif not split_deliverables:
            base.errors.append('At least one deliverable (--deliverables) '
                               'must be specified, but none were entered.')
        if 'bam' in split_deliverables and 'sam' in split_deliverables:
            base.errors.append('Both "bam" and "sam" were entered among '
                               'deliverables (--deliverables), but only '
                               'one should be chosen.')
        base.tsv = 'tsv' in split_deliverables
        base.idx = 'idx' in split_deliverables
        base.bam = 'bam' in split_deliverables or 'sam' in split_deliverables
        base.output_sam = 'sam' in split_deliverables
        base.bed = 'bed' in split_deliverables
        base.bw = 'bw' in split_deliverables
        base.jx = 'jx' in split_deliverables
        base.all_final_outs = (base.tsv or base.bam or base.bed or base.bw)
        base.count_filename = (base.tsv_basename + '.' if base.tsv_basename
                                        else '') + 'counts.tsv.gz'
        if base.bw or base.tsv or base.bam:
            # Outputting average coverages too; put read_counts in output dir
            base.read_counts_out = ab.Url(
                path_join(elastic, base.output_dir, 'cross_sample_results')
            ).to_url(caps=True)
            base.read_counts_file = ab.Url(
                path_join(elastic,
                    base.output_dir,
                    'cross_sample_results',
                    '{0}#{0}'.format(base.count_filename))
            ).to_native_url()
        else:
            # Place in precoverage intermediate directory
            base.read_counts_out = ab.Url(
                path_join(elastic, base.intermediate_dir,
                            'precoverage' if not elastic else '',
                            'counts_for_normalize')
            ).to_url(caps=True)
            base.read_counts_file = ab.Url(
                                    path_join(elastic,
                                    base.intermediate_dir,
                                    'precoverage'
                                    if not elastic
                                    else '',
                                    'counts_for_normalze',
                                    '{0}#{0}'.format(base.count_filename))
                                ).to_native_url()
        isofrag_basename = _transcript_fragment_idx_basename
        if base.isofrag_idx is not None:
            if not (base.isofrag_idx.endswith('.tar.gz') 
                    or base.isofrag_idx.endswith('.tgz')):
                base.errors.append('Isofrag index (--isofrag-idx) must be '
                                   'a compressed archive whose name ends '
                                   'with either ".tgz" or ".tar.gz", but '
                                   '"{0}" was entered'.format(
                                            base.isofrag_idx
                                        )
                                )
            else:
                isofrag_basename = (
                        os.path.basename(base.isofrag_idx)[:-7]
                        if base.isofrag_idx[-7:] == '.tar.gz'
                        else os.path.basename(base.isofrag_idx)[:-4]
                    )
            if not base.all_final_outs:
                base.errors.append('Isofrag index (--isofrag-idx) was '
                                   'specified, so the specified deliverables '
                                   'cannot be obtained. Do not specify '
                                   'the isofrag index to construct a new '
                                   'index.')
            base.transcript_archive = ab.Url(
                                    '{0}#{1}'.format(
                                            base.isofrag_idx,
                                            isofrag_basename
                                        )
                                ).to_native_url()
        elif base.idx:
            # Place transcript index in output directory
            base.transcript_out = ab.Url(
                path_join(elastic, base.output_dir, 'cross_sample_results')
            ).to_url(caps=True)
            base.transcript_archive = ab.Url(
                path_join(elastic,
                    base.output_dir,
                    'cross_sample_results',
                    '{0}.tar.gz#{0}'.format(isofrag_basename))
            ).to_native_url()
        else:
            # Place transcript index in intermediate directory
            base.transcript_out = ab.Url(
                path_join(elastic, base.intermediate_dir,
                            'cojunction_enum' if not elastic else '',
                            isofrag_basename)
            ).to_url(caps=True)
            base.transcript_archive = ab.Url(
                                    path_join(elastic,
                                    base.intermediate_dir,
                                    'cojunction_enum'
                                    if not elastic
                                    else '',
                                    _transcript_fragment_idx_basename,
                                    '{0}.tar.gz#{0}'.format(isofrag_basename))
                                        ).to_native_url()
        base.transcript_in = '{0}/{0}'.format(isofrag_basename)
        # Correct transcript archive in case it's specified by user

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
                default=exe_paths.bowtie1,
                help=('path to Bowtie 1 executable (def: %s)'
                        % (exe_paths.bowtie1
                            if exe_paths.bowtie1 is not None
                            else 'bowtie'))
            )
            exec_parser.add_argument(
                '--bowtie1-build', type=str, required=False,
                metavar='<exe>',
                default=exe_paths.bowtie1_build,
                help=('path to Bowtie 1 Build executable (def: %s)'
                        % (exe_paths.bowtie1_build
                            if exe_paths.bowtie1_build is not None
                            else 'bowtie-build'))
            )
            required_parser.add_argument(
                '-x', '--bowtie-idx', type=str, required=True,
                metavar='<idx | idx/idx>',
                nargs='+',
                help=('Bowtie 1 and 2 index basenames if they share the '
                      'same path; otherwise, Bowtie 1 index basename '
                      'and Bowtie 2 index basename separated by comma or '
                      'space')
            )
            exec_parser.add_argument(
                '--bowtie2', type=str, required=False,
                metavar='<exe>',
                default=exe_paths.bowtie2,
                help=('path to Bowtie 2 executable (def: %s)'
                        % (exe_paths.bowtie2
                            if exe_paths.bowtie2 is not None
                            else 'bowtie2'))
            )
            exec_parser.add_argument(
                '--bowtie2-build', type=str, required=False,
                metavar='<exe>',
                default=exe_paths.bowtie2_build,
                help=('path to Bowtie 2 Build executable (def: %s)'
                        % (exe_paths.bowtie2_build
                            if exe_paths.bowtie2_build is not None
                            else 'bowtie2-build'))
            )
            exec_parser.add_argument(
                '--samtools', type=str, required=False,
                metavar='<exe>',
                default=exe_paths.samtools,
                help=('path to SAMTools executable (def: %s)'
                        % (exe_paths.samtools
                            if exe_paths.samtools is not None
                            else 'samtools'))
            )
            exec_parser.add_argument(
                '--bedgraphtobigwig', type=str, required=False,
                metavar='<exe>',
                default=exe_paths.bedgraphtobigwig,
                help=('path to BedGraphToBigWig executable '
                      '(def: %s)' 
                        % (exe_paths.bedgraphtobigwig
                            if exe_paths.bedgraphtobigwig is not None
                            else 'bedGraphToBigWig'))
            )
            exec_parser.add_argument(
                '--fastq-dump', type=str, required=False,
                metavar='<exe>',
                default=exe_paths.fastq_dump,
                help=('path to SRA Tools fastq-dump executable '
                      '(def: %s)' 
                        % (exe_paths.fastq_dump
                            if exe_paths.fastq_dump is not None
                            else 'fastq-dump'))
            )
            exec_parser.add_argument(
                '--vdb-config', type=str, required=False,
                metavar='<exe>',
                default=exe_paths.vdb_config,
                help=('path to SRA Tools vdb-config executable '
                      '(def: %s)' 
                        % (exe_paths.vdb_config
                            if exe_paths.vdb_config is not None
                            else 'vdb-config'))
            )
        else:
            required_parser.add_argument(
                '-a', '--assembly', type=str, required=True,
                metavar='<choice | tgz>',
                help=('assembly to use for alignment. <choice> can be in '
                      '{"hg19"}. otherwise, specify path to tar.gz Rail '
                      'archive on S3')
            )
        algo_parser.add_argument(
            '-k', type=int, required=False,
            metavar='<int>',
            default=None,
            help=('report up to <int> alignments per read; takes precedence '
                  'over any -k value specified in --bowtie2-args (def: no k)')
        )
        algo_parser.add_argument(
            '--bowtie2-args', type=str, required=False,
            default='',
            metavar='<str>',
            help=('arguments to pass to Bowtie 2, which is always run in '
                  '"--local" mode (def: Bowtie 2 defaults)')
        )
        algo_parser.add_argument(
            '--isofrag-idx', '-fx', type=str, required=False,
            metavar='<file>',
            default=None,
            help=('tar.gz containing transcript fragment (isofrag) index from '
                  'previous Rail-RNA run to use in lieu of searching for '
                  'junctions (def: none; search for junctions)')
        )
        algo_parser.add_argument(
            '--partition-length', type=int, required=False,
            metavar='<int>',
            default=5000,
            help=('smallest unit of genome addressable by single task when '
                  'computing coverage')
        )
        algo_parser.add_argument(
            '--max-readlet-size', type=int, required=False,
            metavar='<int>',
            default=25,
            help=('max size of read segment to align when searching for '
                  'junctions')
        )
        algo_parser.add_argument(
            '--readlet-config-size', type=int, required=False,
            metavar='<int>',
            default=35,
            help=('max number of exonic bases spanned by a path enumerated in '
                  'junction DAG')
        )
        algo_parser.add_argument(
            '--min-readlet-size', type=int, required=False,
            metavar='<int>',
            default=15,
            help=('min size of read segment to align when searching for '
                  'junctions')
        )
        algo_parser.add_argument(
            '--readlet-interval', type=int, required=False,
            metavar='<int>',
            default=4,
            help=('distance between start positions of successive overlapping '
                  'read segments to align when searching for junctions')
        )
        algo_parser.add_argument(
            '--cap-size-multiplier', type=float, required=False,
            default=1.1,
            help=argparse.SUPPRESS
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
            '--do-not-drop-polyA-tails', action='store_const',
            const=True, default=False,
            help=('do not ignore any readlet that is a string of only A '
                  'nucleotides')
        )
        algo_parser.add_argument(
            '--search-filter', type=str, required=False,
            metavar='<choice/int>',
            default='none',
            help=('filter out reads searched for junctions that fall below '
                  'threshold <int> for initially detected anchor length; '
                  'or select <choice> from {"strict", "mild", "none"}')
        )
        algo_parser.add_argument(
            '--junction-criteria', type=str, required=False,
            metavar='<dec,int>',
            default='0.05,5',
            nargs='+',
            help=('if parameter is "f,c", filter out junctions that are not '
                  'either present in at least a fraction f of samples or '
                  'detected in at least c reads of one sample')
        )
        algo_parser.add_argument(
            '--motif-search-window-size', type=int, required=False,
            default=1000,
            help=argparse.SUPPRESS
        )
        algo_parser.add_argument(
            '--max-gaps-mismatches', type=int, required=False,
            default=None,
            help=argparse.SUPPRESS
        )
        algo_parser.add_argument(
            '--motif-radius', type=int, required=False,
            default=5,
            help=argparse.SUPPRESS
        )
        algo_parser.add_argument(
            '--genome-bowtie1-args', type=str, required=False,
            default='-v 0 -a -m 30',
            help=argparse.SUPPRESS
        )
        algo_parser.add_argument(
            '--transcriptome-bowtie2-args', type=str, required=False,
            default='-k 30',
            help=argparse.SUPPRESS
        )
        algo_parser.add_argument(
            '--count-multiplier', type=int, required=False,
            default=15,
            help=argparse.SUPPRESS
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
            '--transcriptome-indexes-per-sample', type=int,
            metavar='<int>',
            default=100,
            help=argparse.SUPPRESS
        )
        algo_parser.add_argument(
            '--library-size', type=int, required=False,
            metavar='<int>',
            default=40,
            help=('Library size (in millions of reads) to which every '
                  'sample\'s coverage is normalized when computing '
                  'average coverage')
        )
        output_parser.add_argument('-d', '--deliverables', required=False,
            metavar='<choice,...>',
            default='idx,tsv,bed,bw',
            nargs='+',
            help=('comma- or space-separated list of desired outputs. Choose '
                  'from among {"idx", "tsv", "bed", "sam" | "bam", "bw", '
                  '"jx"}.')
        )
        output_parser.add_argument(
            '--drop-deletions', action='store_const', const=True,
            default=False,
            help='drop deletions from coverage vectors encoded in bigWigs'
        )
        output_parser.add_argument(
            '--do-not-output-bam-by-chr', action='store_const', const=True,
            default=False,
            help=('place all of a sample\'s alignments in one file rather '
                  'than dividing them up by chromosome')
        )
        output_parser.add_argument(
            '--do-not-output-ave-bw-by-chr', action='store_const',
            const=True, default=False,
            help=('place all of a sample\'s average coverage values in one '
                  'file rather than dividing them up by chromosome')
        )
        output_parser.add_argument(
            '--indel-criteria', type=str, required=False,
            metavar='<dec,int>',
            default='0.05,5',
            nargs='+',
            help=('if parameter is "f,c", suppress indels from cross-sample '
                  'TSVs that are not either present in at least a fraction f '
                  'of samples or detected in at least c reads of one sample')
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
        output_parser.add_argument(
            '--tsv-basename', type=str, required=False,
            metavar='<str>',
            default='',
            help='basename for TSV output (def: *empty*)'
        )
        output_parser.add_argument(
            '--idx-basename', type=str, required=False,
            metavar='<str>',
            default='',
            help=('basename for transcript fragment index output (def: %s)'
                    % _transcript_fragment_idx_basename)
        )

    @staticmethod
    def protosteps(base, input_dir, elastic=False):
        manifest = ('/mnt/space/MANIFEST' if elastic else base.manifest)
        verbose = ('--verbose' if base.verbose else '')
        drop_deletions = ('--drop-deletions' if base.drop_deletions else '')
        keep_alive = ('--keep-alive' if elastic else '')
        scratch  = (('--scratch %s' % base.scratch)
                    if (((not hasattr(base, 'scratch'))
                            or base.scratch is not None)
                         and not elastic) else '')
        output_by_chr = ('--output-bam-by-chr'
                            if not base.do_not_output_bam_by_chr
                            else '')
        realign = (base.bam or base.tsv or base.bed or base.bw)
        nodemanager_mem = (base.nodemanager_mem if hasattr(
                                                        base, 
                                                        'nodemanager_mem'
                                                    )
                            else 1)
        max_tasks = (base.max_tasks if hasattr(
                                                        base, 
                                                        'nodemanager_mem'
                                                    )
                            else 1)
        steps_to_return = [
            {
                'name' : 'Align reads %s' % ('and segment them into readlets'
                                                if base.isofrag_idx is None
                                                else 'to genome'),
                'reducer' : (
                         'align_reads.py --bowtie-idx={0} --bowtie2-idx={1} '
                         '--bowtie2-exe={2} '
                         '--exon-differentials --partition-length={3} '
                         '--min-exon-size={4} '
                         '--search-filter={5} '
                         '--manifest={6} '
                         '--max-readlet-size={7} '
                         '--readlet-interval={8} '
                         '--capping-multiplier={9} '
                         '--gzip-level {10} '
                         '--index-count {11} '
                         '--tie-margin {12} '
                         '{13} {14} {15} {16} {17} {18} {19} -- {20}'
                        ).format(
                                    base.bowtie1_idx,
                                    base.bowtie2_idx,
                                    base.bowtie2_exe,
                                    base.partition_length,
                                    base.min_exon_size,
                                    base.search_filter,
                                    manifest,
                                    base.max_readlet_size,
                                    base.readlet_interval,
                                    base.cap_size_multiplier,
                                    base.gzip_level
                                    if 'gzip_level' in
                                    dir(base) else 3,
                                    base.transcriptome_indexes_per_sample *
                                    base.sample_count,
                                    base.tie_margin,
                                    drop_deletions,
                                    verbose,
                                    keep_alive,
                                    scratch,
                                    output_by_chr,
                                    '--no-realign' if not realign else '',
                                    '--no-polyA'
                                    if not base.do_not_drop_polyA_tails
                                    else '',
                                    base.bowtie2_args + (
                                            ' -p 2 --reorder '
                                            if elastic else ''
                                        ) # 2x threads on EMR cuz reducers/2
                                ),
                'inputs' : [input_dir],
                'no_input_prefix' : True,
                'output' : 'align_reads',
                'tasks' : ('%d,' % (base.sample_count * 3))
                                if elastic else '1x',
                'partition' : '-k1,1',
                'multiple_outputs' : True,
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size * 2),
                        'elephantbird.combined.split.count={task_count}',
                        'mapreduce.reduce.memory.mb=%d'
                        % (nodemanager_mem / max_tasks * 2),
                        'mapreduce.reduce.java.opts=-Xmx%dm'
                        % (nodemanager_mem / max_tasks * 16 / 10)
                    ]
            },
            {
                'name' : 'Align unique readlets to genome',
                'reducer' : (
                         'align_readlets.py --bowtie-idx={0} '
                         '--bowtie-exe={1} {2} {3} --gzip-level={4} {5} '
                         '-- -t --sam-nohead --startverbose {6}').format(
                                                    base.bowtie1_idx,
                                                    base.bowtie1_exe,
                                                    verbose,
                                                    keep_alive,
                                                    base.gzip_level
                                                    if 'gzip_level' in
                                                    dir(base) else 3,
                                                    scratch,
                                                    base.genome_bowtie1_args,
                                                ),
                'inputs' : [path_join(elastic, 'align_reads', 'readletized')],
                'output' : 'align_readlets',
                'tasks' :  ('%d,' % (base.sample_count * 3))
                                if elastic else '1x',
                'partition' : '-k1,1',
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size * 2),
                        'elephantbird.combined.split.count={task_count}'
                    ]
            } if base.isofrag_idx is None else {},
            {
                'name' : 'Search for junctions using readlet alignments',
                'reducer' : (
                         'junction_search.py --bowtie-idx={0} '
                         '--partition-length={1} --max-intron-size={2} '
                         '--min-intron-size={3} --min-exon-size={4} '
                         '--search-window-size={5} {6} '
                         '--motif-radius={7} {8}').format(
                                                base.bowtie1_idx,
                                                base.partition_length,
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
                'output' : 'junction_search',
                'tasks' : ('%d,' % (base.sample_count * 3))
                                    if elastic else '1x',
                'partition' : '-k1,1',
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size * 2),
                        'elephantbird.combined.split.count={task_count}'
                    ]
            } if base.isofrag_idx is None else {},
            {
                'name' : 'Filter out junctions violating confidence criteria',
                'reducer' : (
                         'junction_filter.py --manifest={0} '
                         '--sample-fraction={1} --coverage-threshold={2} '
                         '{3} {4}').format(
                                        manifest,
                                        base.sample_fraction,
                                        base.coverage_threshold,
                                        verbose,
                                        '--collect-junctions'
                                        if base.jx else ''
                                    ),
                'inputs' : ['junction_search'],
                'output' : 'junction_filter',
                'multiple_outputs' : True,
                'tasks' : ('%d,' % max(base.sample_count / 10, 1))
                                    if elastic else '1x',
                'partition' : '-k1,3',
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size),
                        'elephantbird.combined.split.count={task_count}'
                    ]
            } if base.isofrag_idx is None else {},
            {
                'name' : 'Write all detected junctions',
                'reducer' : ('junction_collect.py --out={0} '
                             '--gzip-level {1} {2}').format(
                                                        ab.Url(
                                                            path_join(elastic,
                                                            base.output_dir,
                                                    'collected_junctions')
                                                        ).to_url(caps=True)
                                                        if elastic
                                                        else path_join(elastic,
                                                            base.output_dir,
                                                    'collected_junctions'),
                                                        base.gzip_level
                                                        if 'gzip_level' in
                                                        dir(base) else 3,
                                                        scratch
                                                    ),
                'inputs' : [path_join(elastic, 'junction_filter', 'collect')],
                'output' : 'junction_collect',
                'tasks' : 1,
                'partition' : '-k1,3',
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size),
                        'elephantbird.combined.split.count={task_count}'
                    ]
            } if (base.jx and base.isofrag_idx is None) else {},
            {
                'name' : 'Enumerate junction cooccurrences on readlets',
                'reducer' : ('junction_config.py '
                             '--readlet-size={0} '
                             '--min-overlap-exon-size={1} {2}').format(
                                                    base.readlet_config_size,
                                                    base.min_exon_size,
                                                    verbose
                                                ),
                'inputs' : [path_join(elastic, 'junction_filter', 'filter')],
                'output' : 'junction_config',
                'tasks' : '1x',
                'partition' : '-k1,2',
                'sort' : '-k1,2 -k3,4',
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size),
                        'elephantbird.combined.split.count={task_count}'
                    ]
            } if (base.isofrag_idx is None and (realign or base.idx)) else {},
            {
                'name' : 'Get isofrags for index construction',
                'reducer' : ('junction_fasta.py --bowtie-idx={0} {1}').format(
                                                        base.bowtie1_idx,
                                                        verbose
                                                    ),
                'inputs' : ['junction_config'],
                'output' : 'junction_fasta',
                'tasks' : '1x',
                'partition' : '-k1,3',
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size),
                        'elephantbird.combined.split.count={task_count}'
                    ]
            } if (base.isofrag_idx is None and (realign or base.idx)) else {},
            {
                'name' : 'Build isofrag index',
                'reducer' : ('junction_index.py --bowtie2-build-exe={0} '
                             '--out={1} --basename {2} {3} {4}').format(
                                            base.bowtie2_build_exe,
                                            base.transcript_out,
                                            _transcript_fragment_idx_basename,
                                            keep_alive,
                                            scratch
                                        ),
                'inputs' : ['junction_fasta',
                                path_join(elastic, 'align_reads', 'dummy')],
                'output' : 'junction_index',
                'tasks' : 1,
                'partition' : '-k1,1',
                'sort' : '-k1,1 -k2,2', # ensures ref names in uniform order!
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size),
                        'elephantbird.combined.split.count={task_count}'
                    ]
            } if (base.isofrag_idx is None and (realign or base.idx)) else {},
            {
                'name' : 'Finalize junction cooccurrences on reads',
                'reducer' : (
                         'cojunction_enum.py --bowtie2-idx={0} '
                         '--gzip-level {1} '
                         '--bowtie2-exe={2} {3} {4} --intermediate-dir {5} '
                         '{6} -- {7}').format(
                                            base.transcript_in,
                                            base.gzip_level
                                            if 'gzip_level' in
                                            dir(base) else 3,
                                            base.bowtie2_exe,
                                            verbose,
                                            keep_alive,
                                            ab.Url(
                                                    base.intermediate_dir
                                                ).to_url(caps=True),
                                            scratch,
                                            base.transcriptome_bowtie2_args
                                        ),
                'inputs' : [path_join(elastic, 'align_reads', 'unique')],
                'output' : 'cojunction_enum',
                'tasks' : ('%d,' % (base.sample_count * 10))
                                if elastic else '1x',
                'archives' : base.transcript_archive,
                'partition' : '-k1,1',
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size * 2),
                        'elephantbird.combined.split.count={task_count}'
                    ]
            } if realign else {},
            {
                'name' : 'Get transcriptome elements for read realignment',
                'reducer' : ('cojunction_fasta.py --bowtie-idx={0} '
                             '--index-count {1} {2}').format(
                                                        base.bowtie1_idx,
                                        base.transcriptome_indexes_per_sample *
                                            base.sample_count,
                                                        verbose
                                                    ),
                'inputs' : ['cojunction_enum'],
                'output' : 'cojunction_fasta',
                'tasks' : ('%d,' % (base.sample_count * 10))
                                if elastic else '1x',
                'partition' : '-k1,3',
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size),
                        'elephantbird.combined.split.count={task_count}'
                    ]
            } if realign else {},
            {
                'name' : 'Align reads to transcriptome elements',
                'reducer' : ('realign_reads.py --bowtie2-exe={0} '
                             '--bowtie2-build-exe={1} '
                             '--gzip-level {2} --count-multiplier {3} '
                             '--tie-margin {4} {5} {6} {7} -- {8}').format(
                                            base.bowtie2_exe,
                                            base.bowtie2_build_exe,
                                            base.gzip_level
                                            if 'gzip_level' in
                                            dir(base) else 3,
                                            base.count_multiplier,
                                            base.tie_margin,
                                            verbose,
                                            keep_alive,
                                            scratch,
                                            base.bowtie2_args
                                        ),
                'inputs' : [path_join(elastic, 'align_reads', 'unmapped'),
                            'cojunction_fasta'],
                'mod_partitioner' : True,
                'output' : 'realign_reads',
                # Ensure that a single reducer isn't assigned too much fasta
                'tasks' : ('%d,' % (base.sample_count * 12))
                                 if elastic else '1x',
                'partition' : '-k1,1',
                'sort' : '-k1,1 -k2,3',
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size * 2),
                        'elephantbird.combined.split.count={task_count}'
                    ]
            } if realign else {},
            {
                'name' : 'Collect and compare read alignments',
                'reducer' : ('compare_alignments.py --bowtie-idx={0} '
                             '--partition-length={1} --exon-differentials '
                             '--tie-margin {2} --manifest={3} '
                             '{4} {5} {6} -- {7}').format(
                                            base.bowtie1_idx,
                                            base.partition_length,
                                            base.tie_margin,
                                            manifest,
                                            drop_deletions,
                                            verbose,
                                            output_by_chr,
                                            base.bowtie2_args
                                        ),
                'inputs' : [path_join(elastic, 'align_reads', 'postponed_sam'),
                            'realign_reads'],
                'output' : 'compare_alignments',
                'tasks' : ('%d,' % (base.sample_count * 12))
                                if elastic else '1x',
                'partition' : '-k1,1',
                'multiple_outputs' : True,
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size * 2),
                        'elephantbird.combined.split.count={task_count}',
                        'mapreduce.reduce.memory.mb=%d'
                        % (nodemanager_mem / max_tasks * 2),
                        'mapreduce.reduce.java.opts=-Xmx%dm'
                        % (nodemanager_mem / max_tasks * 16 / 10)
                    ]
            } if realign else {},
            {
                'name' : 'Associate spliced reads with junction coverages',
                'reducer' : 'junction_coverage.py --bowtie-idx {0}'.format(
                                                        base.bowtie1_idx
                                                    ),
                'inputs' : [path_join(elastic, 'compare_alignments',
                                                    'junction_bed'),
                            path_join(elastic, 'compare_alignments',
                                               'sam_junction_ties')],
                'output' : 'junction_coverage',
                'tasks' : '1x',
                'partition' : '-k1,6',
                'sort' : '-k1,6 -k7,7',
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size),
                        'elephantbird.combined.split.count={task_count}'
                    ]
            } if realign else {},
            {
                'name' : 'Finalize primary alignments of spliced reads',
                'reducer' : ('break_ties.py --exon-differentials '
                            '--bowtie-idx {0} --partition-length {1} '
                            '--manifest {2} --tie-margin {3} {4} '
                            '{5} -- {6}').format(
                                    base.bowtie1_idx,
                                    base.partition_length,
                                    manifest,
                                    base.tie_margin,
                                    drop_deletions,
                                    output_by_chr,
                                    base.bowtie2_args
                                ),
                'inputs' : ['junction_coverage',
                            path_join(elastic, 'compare_alignments',
                                               'sam_clip_ties')],
                'output' : 'break_ties',
                'tasks' : '1x',
                'partition' : '-k1,1',
                'multiple_outputs' : True,
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size),
                        'elephantbird.combined.split.count={task_count}'
                    ]
            } if realign else {},
            {
                'name' : (('Write %s with alignments by sample'
                            % ('SAMs' if base.output_sam else 'BAMs'))
                            if base.bam
                            else 'Count mapped reads by contig/sample'),
                'reducer' : (
                         'bam.py --out={0} --bowtie-idx={1} '
                         '--samtools-exe={2} --bam-basename={3} '
                         '--manifest={4} {5} {6} {7} {8} {9} '
                         '--tie-margin {10}').format(
                                        ab.Url(
                                            path_join(elastic,
                                            base.output_dir, 'alignments')
                                        ).to_url(caps=True)
                                        if elastic
                                        else path_join(elastic,
                                            base.output_dir, 'alignments'),
                                        base.bowtie1_idx,
                                        base.samtools_exe,
                                        base.bam_basename,
                                        manifest,
                                        keep_alive,
                                        '--output-by-chromosome'
                                        if (not base.do_not_output_bam_by_chr
                                            or not base.bam) else '',
                                        scratch,
                                        '--output-sam' if base.output_sam
                                        else '',
                                        '--suppress-bam' if not base.bam
                                        else '',
                                        base.tie_margin
                                    ),
                'inputs' : [path_join(elastic, 'compare_alignments', 'sam'),
                            path_join(elastic, 'break_ties', 'sam')]
                            + ([path_join(elastic, 'align_reads', 'sam')]
                                if base.k in [1, None] else []),
                'multiple_outputs' : True,
                'mod_partitioner' : True,
                'output' : 'bam',
                'tasks' : '1x',
                'partition' : '-k1,1',
                'sort' : '-k1,1 -k2,3',
                'extra_args' : [
                        'mapreduce.reduce.shuffle.input.buffer.percent=0.4',
                        'mapreduce.reduce.shuffle.merge.percent=0.4',
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size),
                        'elephantbird.combined.split.count={task_count}'
                    ]
            } if (base.bw or base.tsv or base.bam) else {},
            {
                'name' : 'Write mapped read counts',
                'reducer' : (
                         'collect_read_stats.py --bowtie-idx={0} --out={1} '
                         '--manifest={2} --gzip-level={3} '
                         '--tsv-basename={4} {5} {6}').format(
                                                    base.bowtie1_idx,
                                                    base.read_counts_out,
                                                    manifest,
                                                    base.gzip_level
                                                    if 'gzip_level' in
                                                    dir(base) else 3,
                                                    base.tsv_basename,
                                                    scratch,
                                                    keep_alive
                                                ),
                'inputs' : [path_join(elastic, 'bam', 'counts')],
                'output' : 'read_counts',
                'tasks' : 1,
                'partition' : '-k1,1',
                'sort' : '-k1,1 -k2,2n',
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size),
                        'elephantbird.combined.split.count={task_count}'
                    ]
            } if (base.bw or base.tsv or base.bam) else {},
            {
                'name' : 'Merge exon differentials at same genomic positions',
                'reducer' : 'sum.py {0}'.format(
                                        keep_alive
                                    ),
                'inputs' : [path_join(elastic, 'align_reads', 'exon_diff'),
                            path_join(elastic, 'compare_alignments',
                                               'exon_diff'),
                            path_join(elastic, 'break_ties', 'exon_diff')],
                'output' : 'collapse',
                'tasks' : ('%d,' % (base.sample_count * 12))
                                if elastic else '1x',
                'partition' : '-k1,4',
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size),
                        'elephantbird.combined.split.count={task_count}'
                    ]
            } if base.bw else {},
            {
                'name' : 'Compile sample coverages from exon differentials',
                'reducer' : ('coverage_pre.py --bowtie-idx={0} '
                             '--library-size {1} --read-counts {2} '
                             '--partition-stats --manifest={3} {4}').format(
                                    base.bowtie1_idx,
                                    base.library_size,
                                    base.count_filename,
                                    manifest,
                                    '--output-ave-bigwig-by-chr'
                                    if not base.do_not_output_ave_bw_by_chr
                                    else ''),
                'inputs' : ['collapse'],
                'output' : 'precoverage',
                'tasks' : ('%d,' % (base.sample_count * 12))
                                if elastic else '1x',
                'partition' : '-k1,1',
                'sort' : '-k1,1 -k2,3',
                'files' : base.read_counts_file,
                'multiple_outputs' : True,
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size),
                        'elephantbird.combined.split.count={task_count}'
                    ]
            } if base.bw else {},
            {
                'name' : 'Write bigWigs with exome coverage by sample',
                'reducer' : (
                         'coverage.py --bowtie-idx={0} --percentile={1} '
                         '--out={2} --bigwig-exe={3} '
                         '--manifest={4} {5} {6}').format(base.bowtie1_idx,
                                                     base.normalize_percentile,
                                                     ab.Url(
                                                        path_join(elastic,
                                                        base.output_dir,
                                                        'coverage_bigwigs')
                                                     ).to_url(caps=True)
                                                     if elastic
                                                     else path_join(elastic,
                                                        base.output_dir,
                                                        'coverage_bigwigs'),
                                                     base.bedgraphtobigwig_exe,
                                                     manifest,
                                                     verbose,
                                                     scratch),
                'inputs' : [path_join(elastic, 'precoverage', 'coverage')],
                'output' : 'coverage',
                'mod_partitioner' : True,
                'tasks' : '1x',
                'partition' : '-k1,1',
                'sort' : '-k1,1 -k2,3',
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size),
                        'elephantbird.combined.split.count={task_count}'
                    ]
            } if base.bw else {},
            {
                'name' : 'Aggregate junctions/indels',
                'reducer' : ('bed_pre.py --manifest={0} '
                             '--sample-fraction={1} --coverage-threshold={2} '
                             '{3} {4}').format(
                                        manifest,
                                        base.indel_sample_fraction,
                                        base.indel_coverage_threshold,
                                        verbose,
                                        keep_alive
                                    ),
                'inputs' : [path_join(elastic, 'compare_alignments',
                                               'indel_bed'),
                            path_join(elastic, 'break_ties', 'indel_bed'),
                            path_join(elastic, 'compare_alignments',
                                               'junction_bed'),
                            path_join(elastic, 'break_ties', 'junction_bed')],
                'output' : 'prebed',
                'multiple_outputs' : True,
                'tasks' : ('%d,' % (base.sample_count * 12))
                            if elastic else '1x',
                'partition' : '-k1,5',
                'sort' : '-k1,5 -k6,6n',
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size),
                        'elephantbird.combined.split.count={task_count}'
                    ]
            } if (base.tsv or base.bed) else {},
            {
                'name' : ('Write normalization factors/junctions/indels'
                            if base.tsv else 'Write normalization factors'),
                'reducer' : ('tsv.py --bowtie-idx={0} --out={1} '
                             '--manifest={2} --gzip-level={3} '
                             '--tsv-basename={4} {5} {6}').format(
                                                    base.bowtie1_idx,
                                                    ab.Url(
                                                        path_join(elastic,
                                                        base.output_dir,
                                                    'cross_sample_results')
                                                     ).to_url(caps=True)
                                                    if elastic
                                                    else path_join(elastic,
                                                        base.output_dir,
                                                    'cross_sample_results'),
                                                    manifest,
                                                    base.gzip_level
                                                    if 'gzip_level' in
                                                    dir(base) else 3,
                                                    base.tsv_basename,
                                                    scratch,
                                                    keep_alive
                                                ),
                'inputs' : ['coverage']
                            + ([path_join(elastic, 'prebed', 'collect')]
                                if base.tsv else []),
                'output' : 'tsv',
                'mod_partitioner' : True,
                'tasks' : '1x',
                'partition' : '-k1,1',
                'sort' : '-k1,1 -k2,5',
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size),
                        'elephantbird.combined.split.count={task_count}'
                    ]
            } if (base.tsv or base.bw) else {},
            {
                'name' : 'Write BEDs with junctions/indels by sample',
                'reducer' : (
                         'bed.py --bowtie-idx={0} --out={1} '
                         '--manifest={2} --bed-basename={3} {4} {5}').format(
                                                        base.bowtie1_idx,
                                                        ab.Url(
                                                            path_join(elastic,
                                                            base.output_dir,
                                                        'junctions_and_indels')
                                                         ).to_url(caps=True)
                                                        if elastic
                                                        else path_join(elastic,
                                                            base.output_dir,
                                                        'junctions_and_indels'
                                                        ),
                                                        manifest,
                                                        base.bed_basename,
                                                        scratch,
                                                        keep_alive
                                                    ),
                'inputs' : [path_join(elastic, 'prebed', 'bed')],
                'output' : 'bed',
                'tasks' : '1x',
                'partition' : '-k1,2',
                'sort' : '-k1,2 -k3,5',
                'extra_args' : [
                        'elephantbird.use.combine.input.format=true',
                        'elephantbird.combine.split.size=%d'
                            % (_base_combine_split_size),
                        'elephantbird.combined.split.count={task_count}'
                    ]
            } if base.bed else {}]
        return [step for step in steps_to_return if step != {}]

    @staticmethod
    def bootstrap(base):
        return [
            {
                'Name' : 'Install Rail-RNA and create JAR dependencies',
                'ScriptBootstrapAction' : {
                    'Args' : [
                        base.elastic_rail_path, _jar_target, '-y', '-s'
                    ], # always say yes, only prep dependencies, and symlink  
                    'Path' : base.install_rail_bootstrap
                }
            },
            {
                'Name' : 'Transfer Bowtie indexes to nodes',
                'ScriptBootstrapAction' : {
                    'Args' : [
                        '/mnt/space',
                        base.index_archive
                    ],
                    'Path' : base.copy_index_bootstrap
                }
            },
            {
                'Name' : 'Transfer manifest file to nodes',
                'ScriptBootstrapAction' : {
                    'Args' : [
                        base.manifest,
                        '/mnt/space',
                        'MANIFEST'
                    ],
                    'Path' : base.copy_bootstrap
                }
            }
        ] + ([{
                'Name' : 'Transfer dbGaP key file to nodes',
                'ScriptBootstrapAction' : {
                    'Args' : [
                        base.dbgap_s3_path,
                        '/mnt/space',
                        'DBGAP.ngc'
                    ],
                    'Path' : base.copy_bootstrap
                }
            }] if hasattr(base, 'dbgap_s3_path') else [])

class RailRnaLocalPreprocessJson(object):
    """ Constructs JSON for local mode + preprocess job flow. """
    def __init__(self, manifest, output_dir, isofrag_idx=None,
        intermediate_dir='./intermediate', force=False, aws_exe=None,
        profile='default', region=None, verbose=False,
        nucleotides_per_input=8000000, gzip_input=True,
        do_not_bin_quals=False, short_read_names=False, skip_bad_records=False,
        num_processes=1, gzip_intermediates=False, gzip_level=3,
        sort_memory_cap=(300*1024), max_task_attempts=4, 
        keep_intermediates=False, check_manifest=True,
        scratch=None, sort_exe=None, dbgap_key=None,
        fastq_dump_exe=None, vdb_config_exe=None):
        base = RailRnaErrors(manifest, output_dir, isofrag_idx=isofrag_idx,
            intermediate_dir=intermediate_dir,
            force=force, aws_exe=aws_exe, profile=profile,
            region=region, verbose=verbose,
            max_task_attempts=max_task_attempts, dbgap_key=dbgap_key)
        RailRnaLocal(base, check_manifest=check_manifest,
            num_processes=num_processes, gzip_intermediates=gzip_intermediates,
            gzip_level=gzip_level, sort_memory_cap=sort_memory_cap,
            keep_intermediates=keep_intermediates, scratch=scratch,
            sort_exe=sort_exe, fastq_dump_exe=fastq_dump_exe,
            vdb_config_exe=vdb_config_exe)
        RailRnaPreprocess(base, nucleotides_per_input=nucleotides_per_input,
            gzip_input=gzip_input, do_not_bin_quals=do_not_bin_quals,
            short_read_names=short_read_names,
            skip_bad_records=skip_bad_records)
        raise_runtime_error(base)
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
    def __init__(self, manifest, output_dir, isofrag_idx=None,
        intermediate_dir='./intermediate', force=False, aws_exe=None,
        profile='default', region=None, verbose=False,
        nucleotides_per_input=8000000, gzip_input=True,
        do_not_bin_quals=False, short_read_names=False, skip_bad_records=False,
        num_processes=1, gzip_intermediates=False, gzip_level=3,
        sort_memory_cap=(300*1024), max_task_attempts=4, ipython_profile=None,
        ipcontroller_json=None, scratch=None, direct_write=False,
        keep_intermediates=False, check_manifest=True, sort_exe=None,
        dbgap_key=None, fastq_dump_exe=None, vdb_config_exe=None):
        rc = ipython_client(ipython_profile=ipython_profile,
                                ipcontroller_json=ipcontroller_json)
        base = RailRnaErrors(manifest, output_dir, isofrag_idx=isofrag_idx,
            intermediate_dir=intermediate_dir,
            force=force, aws_exe=aws_exe, profile=profile,
            region=region, verbose=verbose,
            max_task_attempts=max_task_attempts, dbgap_key=None)
        RailRnaLocal(base, check_manifest=check_manifest,
            num_processes=len(rc), gzip_intermediates=gzip_intermediates,
            gzip_level=gzip_level, sort_memory_cap=sort_memory_cap,
            keep_intermediates=keep_intermediates, scratch=scratch,
            direct_write=direct_write, local=False, parallel=False,
            sort_exe=sort_exe, fastq_dump_exe=fastq_dump_exe,
            vdb_config_exe=vdb_config_exe)
        if ab.Url(base.output_dir).is_local:
            '''Add NFS prefix to ensure tasks first copy files to temp dir and
            subsequently upload to final destination.'''
            base.output_dir = ''.join(['nfs://', os.path.abspath(
                                                        base.output_dir
                                                    )])
        RailRnaPreprocess(
                base,
                nucleotides_per_input=nucleotides_per_input,
                gzip_input=gzip_input,
                do_not_bin_quals=do_not_bin_quals,
                short_read_names=short_read_names,
                skip_bad_recoreds=skip_bad_records
            )
        raise_runtime_error(base)
        temp_base_path = ready_engines(rc, base, prep=True)
        engine_bases = {}
        for i in rc.ids:
            engine_bases[i] = RailRnaErrors(
                    manifest, output_dir, isofrag_idx=isofrag_idx,
                    intermediate_dir=intermediate_dir,
                    force=force, aws_exe=aws_exe, profile=profile,
                    region=region, verbose=verbose,
                    max_task_attempts=max_task_attempts,
                    dbgap_key=dbgap_key
                )
        apply_async_with_errors(rc, rc.ids, RailRnaLocal, engine_bases,
            check_manifest=check_manifest, num_processes=num_processes,
            gzip_intermediates=gzip_intermediates, gzip_level=gzip_level,
            sort_memory_cap=sort_memory_cap, scratch=scratch,
            direct_write=direct_write, keep_intermediates=keep_intermediates,
            local=False, parallel=True, ansible=ab.Ansible(),
            sort_exe=sort_exe, fastq_dump_exe=fastq_dump_exe,
            vdb_config_exe=vdb_config_exe)
        engine_base_checks = {}
        for i in rc.ids:
            engine_base_checks[i] = engine_bases[i].check_program
        if base.check_curl_on_engines:
            apply_async_with_errors(rc, rc.ids, engine_base_checks,
                'curl', 'cURL', '--curl', entered_exe=base.curl_exe,
                reason=base.check_curl_on_engines, is_exe=is_exe, which=which)
        engine_base_checks = {}
        for i in rc.ids:
            engine_base_checks[i] = engine_bases[i].check_s3
        if base.check_s3_on_engines:
            apply_async_with_errors(rc, rc.ids, engine_base_checks,
                reason=base.check_curl_on_engines, is_exe=is_exe, which=which)
        raise_runtime_error(base)
        self._json_serial = {}
        step_dir = os.path.join(temp_base_path, 'rna', 'steps')
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
    def __init__(self, manifest, output_dir, isofrag_idx=None,
        intermediate_dir='./intermediate', force=False, aws_exe=None,
        profile='default', region=None,
        service_role=None, instance_profile=None,
        verbose=False, nucleotides_per_input=8000000, gzip_input=True,
        do_not_bin_quals=False, short_read_names=False, skip_bad_records=False,
        log_uri=None, ami_version='3.8.0',
        visible_to_all_users=False, tags='',
        name='Rail-RNA Job Flow',
        action_on_failure='TERMINATE_JOB_FLOW',
        hadoop_jar=None,
        master_instance_count=1, master_instance_type='c1.xlarge',
        master_instance_bid_price=None, core_instance_count=1,
        core_instance_type=None, core_instance_bid_price=None,
        task_instance_count=0, task_instance_type=None,
        task_instance_bid_price=None, ec2_key_name=None,
        ec2_subnet_id=None, ec2_master_security_group_id=None,
        ec2_slave_security_group_id=None, keep_alive=False,
        termination_protected=False, consistent_view=False,
        no_direct_copy=False, check_manifest=True, intermediate_lifetime=4,
        max_task_attempts=4, secure_stack_name=None, dbgap_key=None):
        base = RailRnaErrors(manifest, output_dir, isofrag_idx=isofrag_idx,
            intermediate_dir=intermediate_dir,
            force=force, aws_exe=aws_exe, profile=profile,
            region=region, service_role=service_role,
            instance_profile=instance_profile, verbose=verbose,
            max_task_attempts=max_task_attempts, dbgap_key=dbgap_key)
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
            ec2_key_name=ec2_key_name, ec2_subnet_id=ec2_subnet_id,
            ec2_master_security_group_id=ec2_master_security_group_id,
            ec2_slave_security_group_id=ec2_slave_security_group_id,
            keep_alive=keep_alive,
            termination_protected=termination_protected,
            consistent_view=consistent_view,
            no_direct_copy=no_direct_copy,
            intermediate_lifetime=intermediate_lifetime,
            secure_stack_name=secure_stack_name)
        RailRnaPreprocess(base, nucleotides_per_input=nucleotides_per_input,
            gzip_input=gzip_input, do_not_bin_quals=do_not_bin_quals,
            short_read_names=short_read_names,
            skip_bad_records=skip_bad_records)
        raise_runtime_error(base)
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
                    base.hadoop_jar, _elastic_step_dir,
                    reducer_count, base.intermediate_dir, unix=True,
                    no_direct_copy=base.no_direct_copy
                )
        self._json_serial['AmiVersion'] = base.ami_version
        self._json_serial['ServiceRole'] = base.service_role
        self._json_serial['JobFlowRole'] = base.instance_profile
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
        self._json_serial['BootstrapActions'] = (
                RailRnaElastic.prebootstrap(base)
                + RailRnaPreprocess.bootstrap(base)
                + RailRnaPreprocess.srabootstrap(base)
                + RailRnaElastic.bootstrap(base)
            )
        self.base = base
    
    @property
    def json_serial(self):
        return self._json_serial

class RailRnaLocalAlignJson(object):
    """ Constructs JSON for local mode + align job flow. """
    def __init__(self, manifest, output_dir, input_dir,
        isofrag_idx=None, intermediate_dir='./intermediate',
        force=False, aws_exe=None, profile='default', region=None,
        verbose=False, bowtie1_exe=None,
        bowtie_idx='genome', bowtie1_build_exe=None, bowtie2_exe=None,
        bowtie2_build_exe=None, k=1, bowtie2_args='',
        samtools_exe=None, bedgraphtobigwig_exe=None,
        partition_length=5000, max_readlet_size=25,
        readlet_config_size=32, min_readlet_size=15, readlet_interval=4,
        cap_size_multiplier=1.2, max_intron_size=500000, min_intron_size=10,
        min_exon_size=9, library_size=40, search_filter='none',
        motif_search_window_size=1000, max_gaps_mismatches=3,
        motif_radius=5, genome_bowtie1_args='-v 0 -a -m 80',
        transcriptome_bowtie2_args='-k 30', count_multiplier=15,
        junction_criteria='0.5,5', indel_criteria='0.5,5', tie_margin=6,
        transcriptome_indexes_per_sample=500, normalize_percentile=0.75,
        drop_deletions=False, do_not_output_bam_by_chr=False,
        do_not_output_ave_bw_by_chr=False, do_not_drop_polyA_tails=False,
        deliverables='idx,tsv,bed,bw', bam_basename='alignments',
        bed_basename='', tsv_basename='', num_processes=1,
        gzip_intermediates=False, gzip_level=3, sort_memory_cap=(300*1024),
        max_task_attempts=4, keep_intermediates=False, scratch=None,
        sort_exe=None):
        base = RailRnaErrors(manifest, output_dir, isofrag_idx=isofrag_idx,
            intermediate_dir=intermediate_dir,
            force=force, aws_exe=aws_exe, profile=profile,
            region=region, verbose=verbose,
            max_task_attempts=max_task_attempts)
        RailRnaLocal(base, check_manifest=False, num_processes=num_processes,
            gzip_intermediates=gzip_intermediates, gzip_level=gzip_level,
            sort_memory_cap=sort_memory_cap,
            keep_intermediates=keep_intermediates,
            scratch=scratch, sort_exe=sort_exe)
        RailRnaAlign(base, input_dir=input_dir,
            elastic=False, bowtie1_exe=bowtie1_exe,
            bowtie_idx=bowtie_idx, bowtie1_build_exe=bowtie1_build_exe,
            bowtie2_exe=bowtie2_exe, bowtie2_build_exe=bowtie2_build_exe,
            k=k, bowtie2_args=bowtie2_args, samtools_exe=samtools_exe,
            bedgraphtobigwig_exe=bedgraphtobigwig_exe,
            partition_length=partition_length,
            max_readlet_size=max_readlet_size,
            readlet_config_size=readlet_config_size,
            min_readlet_size=min_readlet_size,
            readlet_interval=readlet_interval,
            cap_size_multiplier=cap_size_multiplier,
            max_intron_size=max_intron_size,
            min_intron_size=min_intron_size, min_exon_size=min_exon_size,
            library_size=library_size,
            search_filter=search_filter,
            motif_search_window_size=motif_search_window_size,
            max_gaps_mismatches=max_gaps_mismatches,
            motif_radius=motif_radius,
            genome_bowtie1_args=genome_bowtie1_args,
            junction_criteria=junction_criteria,
            indel_criteria=indel_criteria,
            transcriptome_bowtie2_args=transcriptome_bowtie2_args,
            count_multiplier=count_multiplier,
            tie_margin=tie_margin,
            normalize_percentile=normalize_percentile,
            transcriptome_indexes_per_sample=transcriptome_indexes_per_sample,
            drop_deletions=drop_deletions,
            do_not_output_bam_by_chr=do_not_output_bam_by_chr,
            do_not_output_ave_bw_by_chr=do_not_output_ave_bw_by_chr,
            do_not_drop_polyA_tails=do_not_drop_polyA_tails,
            deliverables=deliverables, bam_basename=bam_basename,
            bed_basename=bed_basename, tsv_basename=tsv_basename)
        raise_runtime_error(base)
        print_to_screen(base.detect_message)
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
        isofrag_idx=None, intermediate_dir='./intermediate',
        force=False, aws_exe=None, profile='default', region=None,
        verbose=False, bowtie1_exe=None,
        bowtie_idx='genome', bowtie1_build_exe=None, bowtie2_exe=None,
        bowtie2_build_exe=None, k=1, bowtie2_args='',
        samtools_exe=None, bedgraphtobigwig_exe=None,
        partition_length=5000, max_readlet_size=25,
        readlet_config_size=32, min_readlet_size=15, readlet_interval=4,
        cap_size_multiplier=1.2, max_intron_size=500000, min_intron_size=10,
        min_exon_size=9, library_size=40, search_filter='none',
        motif_search_window_size=1000, max_gaps_mismatches=3,
        motif_radius=5, genome_bowtie1_args='-v 0 -a -m 80',
        transcriptome_bowtie2_args='-k 30', count_multiplier=15,
        junction_criteria='0.5,5', indel_criteria='0.5,5', tie_margin=6,
        transcriptome_indexes_per_sample=500, normalize_percentile=0.75,
        drop_deletions=False, do_not_output_bam_by_chr=False,
        do_not_output_ave_bw_by_chr=False, do_not_drop_polyA_tails=False,
        deliverables='idx,tsv,bed,bw', bam_basename='alignments',
        bed_basename='', tsv_basename='', num_processes=1,
        ipython_profile=None, ipcontroller_json=None, scratch=None,
        direct_write=False, gzip_intermediates=False, gzip_level=3,
        sort_memory_cap=(300*1024), max_task_attempts=4,
        keep_intermediates=False, do_not_copy_index_to_nodes=False,
        sort_exe=None):
        rc = ipython_client(ipython_profile=ipython_profile,
                                ipcontroller_json=ipcontroller_json)
        base = RailRnaErrors(manifest, output_dir, isofrag_idx=isofrag_idx,
            intermediate_dir=intermediate_dir,
            force=force, aws_exe=aws_exe, profile=profile,
            region=region, verbose=verbose,
            max_task_attempts=max_task_attempts)
        RailRnaLocal(base, check_manifest=False,
            num_processes=len(rc), gzip_intermediates=gzip_intermediates,
            gzip_level=gzip_level, sort_memory_cap=sort_memory_cap,
            keep_intermediates=keep_intermediates,
            local=False, parallel=False, scratch=scratch,
            direct_write=direct_write, sort_exe=sort_exe)
        if ab.Url(base.output_dir).is_local:
            '''Add NFS prefix to ensure tasks first copy files to temp dir and
            subsequently upload to S3.'''
            base.output_dir = ''.join(['nfs://', os.path.abspath(
                                                        base.output_dir
                                                    )])
        RailRnaAlign(base, input_dir=input_dir,
            elastic=False, bowtie1_exe=bowtie1_exe,
            bowtie_idx=bowtie_idx, bowtie1_build_exe=bowtie1_build_exe,
            bowtie2_exe=bowtie2_exe, bowtie2_build_exe=bowtie2_build_exe,
            k=k, bowtie2_args=bowtie2_args, samtools_exe=samtools_exe,
            bedgraphtobigwig_exe=bedgraphtobigwig_exe,
            partition_length=partition_length,
            max_readlet_size=max_readlet_size,
            readlet_config_size=readlet_config_size,
            min_readlet_size=min_readlet_size,
            readlet_interval=readlet_interval,
            cap_size_multiplier=cap_size_multiplier,
            max_intron_size=max_intron_size,
            min_intron_size=min_intron_size, min_exon_size=min_exon_size,
            library_size=library_size,
            search_filter=search_filter,
            motif_search_window_size=motif_search_window_size,
            max_gaps_mismatches=max_gaps_mismatches,
            motif_radius=motif_radius,
            genome_bowtie1_args=genome_bowtie1_args,
            junction_criteria=junction_criteria,
            indel_criteria=indel_criteria,
            transcriptome_bowtie2_args=transcriptome_bowtie2_args,
            count_multiplier=count_multiplier,
            tie_margin=tie_margin,
            normalize_percentile=normalize_percentile,
            transcriptome_indexes_per_sample=transcriptome_indexes_per_sample,
            drop_deletions=drop_deletions,
            do_not_output_bam_by_chr=do_not_output_bam_by_chr,
            do_not_output_ave_bw_by_chr=do_not_output_ave_bw_by_chr,
            do_not_drop_polyA_tails=do_not_drop_polyA_tails,
            deliverables=deliverables, bam_basename=bam_basename,
            tsv_basename=tsv_basename, bed_basename=bed_basename)
        raise_runtime_error(base)
        temp_base_path = ready_engines(rc, base, prep=False)
        engine_bases = {}
        for i in rc.ids:
            engine_bases[i] = RailRnaErrors(
                    manifest, output_dir, isofrag_idx=isofrag_idx,
                    intermediate_dir=intermediate_dir,
                    force=force, aws_exe=aws_exe, profile=profile,
                    region=region, verbose=verbose,
                    max_task_attempts=max_task_attempts
                )
        apply_async_with_errors(rc, rc.ids, RailRnaLocal, engine_bases,
            check_manifest=False, num_processes=num_processes,
            gzip_intermediates=gzip_intermediates, gzip_level=gzip_level,
            sort_memory_cap=sort_memory_cap,
            keep_intermediates=keep_intermediates, local=False, parallel=True,
            ansible=ab.Ansible(), scratch=scratch, direct_write=direct_write,
            sort_exe=sort_exe)
        apply_async_with_errors(rc, rc.ids, RailRnaAlign, engine_bases,
            input_dir=input_dir, elastic=False, bowtie1_exe=bowtie1_exe,
            bowtie_idx=bowtie_idx, bowtie1_build_exe=bowtie1_build_exe,
            bowtie2_exe=bowtie2_exe, bowtie2_build_exe=bowtie2_build_exe,
            k=k, bowtie2_args=bowtie2_args, samtools_exe=samtools_exe,
            bedgraphtobigwig_exe=bedgraphtobigwig_exe,
            partition_length=partition_length,
            max_readlet_size=max_readlet_size,
            readlet_config_size=readlet_config_size,
            min_readlet_size=min_readlet_size,
            readlet_interval=readlet_interval,
            cap_size_multiplier=cap_size_multiplier,
            max_intron_size=max_intron_size,
            min_intron_size=min_intron_size, min_exon_size=min_exon_size,
            library_size=library_size, 
            search_filter=search_filter,
            motif_search_window_size=motif_search_window_size,
            max_gaps_mismatches=max_gaps_mismatches,
            motif_radius=motif_radius,
            genome_bowtie1_args=genome_bowtie1_args,
            junction_criteria=junction_criteria,
            indel_criteria=indel_criteria,
            transcriptome_bowtie2_args=transcriptome_bowtie2_args,
            count_multiplier=count_multiplier,
            tie_margin=tie_margin,
            normalize_percentile=normalize_percentile,
            transcriptome_indexes_per_sample=transcriptome_indexes_per_sample,
            drop_deletions=drop_deletions,
            do_not_output_bam_by_chr=do_not_output_bam_by_chr,
            do_not_output_ave_bw_by_chr=do_not_output_ave_bw_by_chr,
            do_not_drop_polyA_tails=do_not_drop_polyA_tails,
            deliverables=deliverables, bam_basename=bam_basename,
            tsv_basename=tsv_basename, bed_basename=bed_basename)
        engine_base_checks = {}
        for i in rc.ids:
            engine_base_checks[i] = engine_bases[i].check_program
        if base.check_curl_on_engines:
            apply_async_with_errors(rc, rc.ids, engine_base_checks,
                'curl', 'cURL', '--curl', entered_exe=base.curl_exe,
                reason=base.check_curl_on_engines, is_exe=is_exe, which=which)
        engine_base_checks = {}
        for i in rc.ids:
            engine_base_checks[i] = engine_bases[i].check_s3
        if base.check_s3_on_engines:
            apply_async_with_errors(rc, rc.ids, engine_base_checks,
                reason=base.check_curl_on_engines, is_exe=is_exe, which=which)
        raise_runtime_error(engine_bases)
        print_to_screen(base.detect_message)
        self._json_serial = {}
        step_dir = os.path.join(temp_base_path, 'rna', 'steps')
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
        isofrag_idx=None, intermediate_dir='./intermediate',
        force=False, aws_exe=None, profile='default', region=None,
        service_role=None, instance_profile=None,
        verbose=False, bowtie1_exe=None, bowtie_idx='genome',
        bowtie1_build_exe=None, bowtie2_exe=None,
        bowtie2_build_exe=None, k=1, bowtie2_args='',
        samtools_exe=None, bedgraphtobigwig_exe=None,
        partition_length=5000, max_readlet_size=25,
        readlet_config_size=32, min_readlet_size=15, readlet_interval=4,
        cap_size_multiplier=1.2, max_intron_size=500000, min_intron_size=10,
        min_exon_size=9, library_size=40, search_filter='none',
        motif_search_window_size=1000, max_gaps_mismatches=3,
        motif_radius=5, genome_bowtie1_args='-v 0 -a -m 80',
        transcriptome_bowtie2_args='-k 30', count_multiplier=15,
        junction_criteria='0.5,5', indel_criteria='0.5,5', tie_margin=6,
        transcriptome_indexes_per_sample=500, normalize_percentile=0.75,
        drop_deletions=False, do_not_output_bam_by_chr=False,
        do_not_output_ave_bw_by_chr=False, do_not_drop_polyA_tails=False,
        deliverables='idx,tsv,bed,bw', bam_basename='alignments',
        bed_basename='', tsv_basename='', log_uri=None, ami_version='3.8.0',
        visible_to_all_users=False, tags='', name='Rail-RNA Job Flow',
        action_on_failure='TERMINATE_JOB_FLOW', hadoop_jar=None,
        master_instance_count=1, master_instance_type='c1.xlarge',
        master_instance_bid_price=None, core_instance_count=1,
        core_instance_type=None, core_instance_bid_price=None,
        task_instance_count=0, task_instance_type=None,
        task_instance_bid_price=None, ec2_key_name=None,
        ec2_subnet_id=None, ec2_master_security_group_id=None,
        ec2_slave_security_group_id=None, keep_alive=False,
        termination_protected=False, consistent_view=False,
        no_direct_copy=False, intermediate_lifetime=4, max_task_attempts=4,
        secure_stack_name=None):
        base = RailRnaErrors(manifest, output_dir, isofrag_idx=isofrag_idx,
            intermediate_dir=intermediate_dir,
            force=force, aws_exe=aws_exe, profile=profile,
            region=region, service_role=service_role,
            instance_profile=instance_profile, verbose=verbose,
            max_task_attempts=max_task_attempts)
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
            ec2_key_name=ec2_key_name, ec2_subnet_id=ec2_subnet_id,
            ec2_master_security_group_id=ec2_master_security_group_id,
            ec2_slave_security_group_id=ec2_slave_security_group_id,
            keep_alive=keep_alive,
            termination_protected=termination_protected,
            consistent_view=consistent_view,
            no_direct_copy=no_direct_copy,
            intermediate_lifetime=intermediate_lifetime,
            secure_stack_name=secure_stack_name)
        RailRnaAlign(base, input_dir=input_dir,
            elastic=True, bowtie1_exe=bowtie1_exe,
            bowtie_idx=bowtie_idx, bowtie1_build_exe=bowtie1_build_exe,
            bowtie2_exe=bowtie2_exe, bowtie2_build_exe=bowtie2_build_exe,
            k=k, bowtie2_args=bowtie2_args, samtools_exe=samtools_exe,
            bedgraphtobigwig_exe=bedgraphtobigwig_exe,
            partition_length=partition_length,
            max_readlet_size=max_readlet_size,
            readlet_config_size=readlet_config_size,
            min_readlet_size=min_readlet_size,
            readlet_interval=readlet_interval,
            cap_size_multiplier=cap_size_multiplier,
            max_intron_size=max_intron_size,
            min_intron_size=min_intron_size, min_exon_size=min_exon_size,
            library_size=library_size, 
            search_filter=search_filter,
            motif_search_window_size=motif_search_window_size,
            max_gaps_mismatches=max_gaps_mismatches,
            motif_radius=motif_radius,
            genome_bowtie1_args=genome_bowtie1_args,
            junction_criteria=junction_criteria,
            indel_criteria=indel_criteria,
            transcriptome_bowtie2_args=transcriptome_bowtie2_args,
            count_multiplier=count_multiplier,
            tie_margin=tie_margin,
            normalize_percentile=normalize_percentile,
            transcriptome_indexes_per_sample=transcriptome_indexes_per_sample,
            drop_deletions=drop_deletions,
            do_not_output_bam_by_chr=do_not_output_bam_by_chr,
            do_not_output_ave_bw_by_chr=do_not_output_ave_bw_by_chr,
            do_not_drop_polyA_tails=do_not_drop_polyA_tails,
            deliverables=deliverables, bam_basename=bam_basename,
            tsv_basename=tsv_basename, bed_basename=bed_basename,
            s3_ansible=ab.S3Ansible(aws_exe=base.aws_exe,
                                        profile=base.profile))
        raise_runtime_error(base)
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
                    base.hadoop_jar, _elastic_step_dir,
                    reducer_count, base.intermediate_dir, unix=True,
                    no_direct_copy=base.no_direct_copy
                )
        self._json_serial['AmiVersion'] = base.ami_version
        self._json_serial['ServiceRole'] = base.service_role
        self._json_serial['JobFlowRole'] = base.instance_profile
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
        self._json_serial['BootstrapActions'] = (
                RailRnaElastic.prebootstrap(base)
                + RailRnaAlign.bootstrap(base)
                + RailRnaElastic.bootstrap(base)
            )
        self.base = base
    
    @property
    def json_serial(self):
        return self._json_serial

class RailRnaLocalAllJson(object):
    """ Constructs JSON for local mode + preprocess+align job flow. """
    def __init__(self, manifest, output_dir, isofrag_idx=None,
        intermediate_dir='./intermediate', force=False, aws_exe=None,
        profile='default', region=None,
        verbose=False, nucleotides_per_input=8000000, gzip_input=True,
        do_not_bin_quals=False, short_read_names=False, skip_bad_records=False,
        bowtie1_exe=None, bowtie_idx='genome', bowtie1_build_exe=None,
        bowtie2_exe=None, bowtie2_build_exe=None, k=1, bowtie2_args='',
        samtools_exe=None, bedgraphtobigwig_exe=None,
        partition_length=5000, max_readlet_size=25,
        readlet_config_size=32, min_readlet_size=15, readlet_interval=4,
        cap_size_multiplier=1.2, max_intron_size=500000, min_intron_size=10,
        min_exon_size=9, library_size=40, search_filter='none',
        motif_search_window_size=1000, max_gaps_mismatches=3,
        motif_radius=5, genome_bowtie1_args='-v 0 -a -m 80',
        junction_criteria='0.5,5', indel_criteria='0.5,5',
        transcriptome_bowtie2_args='-k 30', tie_margin=6, count_multiplier=15,
        transcriptome_indexes_per_sample=500, normalize_percentile=0.75,
        drop_deletions=False, do_not_output_bam_by_chr=False,
        do_not_output_ave_bw_by_chr=False, do_not_drop_polyA_tails=False,
        deliverables='idx,tsv,bed,bw', bam_basename='alignments',
        bed_basename='', tsv_basename='', num_processes=1,
        gzip_intermediates=False, gzip_level=3,
        sort_memory_cap=(300*1024), max_task_attempts=4,
        keep_intermediates=False, check_manifest=True, scratch=None,
        sort_exe=None, dbgap_key=None, fastq_dump_exe=None,
        vdb_config_exe=None):
        base = RailRnaErrors(manifest, output_dir, isofrag_idx=isofrag_idx,
            intermediate_dir=intermediate_dir,
            force=force, aws_exe=aws_exe, profile=profile,
            region=region, verbose=verbose,
            max_task_attempts=max_task_attempts, dbgap_key=dbgap_key)
        RailRnaPreprocess(base, nucleotides_per_input=nucleotides_per_input,
            gzip_input=gzip_input, do_not_bin_quals=do_not_bin_quals,
            short_read_names=short_read_names,
            skip_bad_records=skip_bad_records)
        RailRnaLocal(base, check_manifest=check_manifest,
            num_processes=num_processes, gzip_intermediates=gzip_intermediates,
            gzip_level=gzip_level, sort_memory_cap=sort_memory_cap,
            keep_intermediates=keep_intermediates, scratch=scratch,
            sort_exe=sort_exe, fastq_dump_exe=fastq_dump_exe,
            vdb_config_exe=vdb_config_exe)
        RailRnaAlign(base, bowtie1_exe=bowtie1_exe,
            bowtie_idx=bowtie_idx, bowtie1_build_exe=bowtie1_build_exe,
            bowtie2_exe=bowtie2_exe, bowtie2_build_exe=bowtie2_build_exe,
            k=k, bowtie2_args=bowtie2_args, samtools_exe=samtools_exe,
            bedgraphtobigwig_exe=bedgraphtobigwig_exe,
            partition_length=partition_length,
            max_readlet_size=max_readlet_size,
            readlet_config_size=readlet_config_size,
            min_readlet_size=min_readlet_size,
            readlet_interval=readlet_interval,
            cap_size_multiplier=cap_size_multiplier,
            max_intron_size=max_intron_size,
            min_intron_size=min_intron_size, min_exon_size=min_exon_size,
            library_size=library_size,
            search_filter=search_filter,
            motif_search_window_size=motif_search_window_size,
            max_gaps_mismatches=max_gaps_mismatches,
            motif_radius=motif_radius,
            genome_bowtie1_args=genome_bowtie1_args,
            transcriptome_bowtie2_args=transcriptome_bowtie2_args,
            count_multiplier=count_multiplier,
            junction_criteria=junction_criteria,
            indel_criteria=indel_criteria,
            tie_margin=tie_margin,
            normalize_percentile=normalize_percentile,
            transcriptome_indexes_per_sample=transcriptome_indexes_per_sample,
            drop_deletions=drop_deletions,
            do_not_output_bam_by_chr=do_not_output_bam_by_chr,
            do_not_output_ave_bw_by_chr=do_not_output_ave_bw_by_chr,
            do_not_drop_polyA_tails=do_not_drop_polyA_tails,
            deliverables=deliverables, bam_basename=bam_basename,
            bed_basename=bed_basename, tsv_basename=tsv_basename)
        raise_runtime_error(base)
        print_to_screen(base.detect_message)
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
    def __init__(self, manifest, output_dir, isofrag_idx=None,
        intermediate_dir='./intermediate', force=False, aws_exe=None,
        profile='default', region=None,
        verbose=False, nucleotides_per_input=8000000, gzip_input=True,
        do_not_bin_quals=False, short_read_names=False, skip_bad_records=False,
        bowtie1_exe=None, bowtie_idx='genome', bowtie1_build_exe=None,
        bowtie2_exe=None, bowtie2_build_exe=None, k=1, bowtie2_args='',
        samtools_exe=None, bedgraphtobigwig_exe=None,
        partition_length=5000, max_readlet_size=25,
        readlet_config_size=32, min_readlet_size=15, readlet_interval=4,
        cap_size_multiplier=1.2, max_intron_size=500000, min_intron_size=10,
        min_exon_size=9, library_size=40, search_filter='none',
        motif_search_window_size=1000, max_gaps_mismatches=3,
        motif_radius=5, genome_bowtie1_args='-v 0 -a -m 80',
        junction_criteria='0.5,5', indel_criteria='0.5,5',
        transcriptome_bowtie2_args='-k 30', tie_margin=6, count_multiplier=15,
        transcriptome_indexes_per_sample=500, normalize_percentile=0.75,
        drop_deletions=False, do_not_output_bam_by_chr=False,
        do_not_output_ave_bw_by_chr=False, do_not_drop_polyA_tails=False,
        deliverables='idx,tsv,bed,bw', bam_basename='alignments',
        bed_basename='', tsv_basename='', num_processes=1,
        gzip_intermediates=False, gzip_level=3, sort_memory_cap=(300*1024),
        max_task_attempts=4, ipython_profile=None, ipcontroller_json=None,
        scratch=None, direct_write=False, keep_intermediates=False,
        check_manifest=True, do_not_copy_index_to_nodes=False, sort_exe=None,
        dbgap_key=None, fastq_dump_exe=None, vdb_config_exe=None):
        rc = ipython_client(ipython_profile=ipython_profile,
                                ipcontroller_json=ipcontroller_json)
        base = RailRnaErrors(manifest, output_dir, isofrag_idx=isofrag_idx,
            intermediate_dir=intermediate_dir,
            force=force, aws_exe=aws_exe, profile=profile,
            region=region, verbose=verbose,
            max_task_attempts=max_task_attempts, dbgap_key=dbgap_key)
        RailRnaLocal(base, check_manifest=check_manifest,
            num_processes=len(rc), gzip_intermediates=gzip_intermediates,
            gzip_level=gzip_level, sort_memory_cap=sort_memory_cap,
            direct_write=direct_write, keep_intermediates=keep_intermediates,
            local=False, parallel=False, scratch=scratch, sort_exe=sort_exe,
            fastq_dump_exe=fastq_dump_exe, vdb_config_exe=vdb_config_exe)
        if ab.Url(base.output_dir).is_local:
            '''Add NFS prefix to ensure tasks first copy files to temp dir and
            subsequently upload to S3.'''
            base.output_dir = ''.join(['nfs://', os.path.abspath(
                                                        base.output_dir
                                                    )])
        RailRnaPreprocess(base, nucleotides_per_input=nucleotides_per_input,
            gzip_input=gzip_input, do_not_bin_quals=do_not_bin_quals,
            short_read_names=short_read_names,
            skip_bad_records=skip_bad_records)
        RailRnaAlign(base, bowtie1_exe=bowtie1_exe,
            bowtie_idx=bowtie_idx, bowtie1_build_exe=bowtie1_build_exe,
            bowtie2_exe=bowtie2_exe, bowtie2_build_exe=bowtie2_build_exe,
            k=k, bowtie2_args=bowtie2_args, samtools_exe=samtools_exe,
            bedgraphtobigwig_exe=bedgraphtobigwig_exe,
            partition_length=partition_length,
            max_readlet_size=max_readlet_size,
            readlet_config_size=readlet_config_size,
            min_readlet_size=min_readlet_size,
            readlet_interval=readlet_interval,
            cap_size_multiplier=cap_size_multiplier,
            max_intron_size=max_intron_size,
            min_intron_size=min_intron_size, min_exon_size=min_exon_size,
            library_size=library_size, 
            search_filter=search_filter,
            motif_search_window_size=motif_search_window_size,
            max_gaps_mismatches=max_gaps_mismatches,
            motif_radius=motif_radius,
            genome_bowtie1_args=genome_bowtie1_args,
            transcriptome_bowtie2_args=transcriptome_bowtie2_args,
            count_multiplier=count_multiplier,
            junction_criteria=junction_criteria,
            indel_criteria=indel_criteria,
            tie_margin=tie_margin,
            normalize_percentile=normalize_percentile,
            transcriptome_indexes_per_sample=transcriptome_indexes_per_sample,
            drop_deletions=drop_deletions,
            do_not_output_bam_by_chr=do_not_output_bam_by_chr,
            do_not_output_ave_bw_by_chr=do_not_output_ave_bw_by_chr,
            do_not_drop_polyA_tails=do_not_drop_polyA_tails,
            deliverables=deliverables, bam_basename=bam_basename,
            bed_basename=bed_basename, tsv_basename=tsv_basename)
        raise_runtime_error(base)
        temp_base_path = ready_engines(rc, base, prep=False)
        engine_bases = {}
        for i in rc.ids:
            engine_bases[i] = RailRnaErrors(
                    manifest, output_dir, isofrag_idx=isofrag_idx,
                    intermediate_dir=intermediate_dir,
                    force=force, aws_exe=aws_exe, profile=profile,
                    region=region, verbose=verbose,
                    max_task_attempts=max_task_attempts,
                    dbgap_key=dbgap_key
                )
        apply_async_with_errors(rc, rc.ids, RailRnaLocal, engine_bases,
            check_manifest=check_manifest, num_processes=num_processes,
            gzip_intermediates=gzip_intermediates, gzip_level=gzip_level,
            sort_memory_cap=sort_memory_cap,
            keep_intermediates=keep_intermediates, local=False, parallel=True,
            ansible=ab.Ansible(), scratch=scratch, direct_write=direct_write,
            sort_exe=sort_exe, fastq_dump_exe=fastq_dump_exe,
            vdb_config_exe=vdb_config_exe)
        apply_async_with_errors(rc, rc.ids, RailRnaAlign, engine_bases,
            bowtie1_exe=bowtie1_exe,
            bowtie_idx=bowtie_idx, bowtie1_build_exe=bowtie1_build_exe,
            bowtie2_exe=bowtie2_exe, bowtie2_build_exe=bowtie2_build_exe,
            k=k, bowtie2_args=bowtie2_args, samtools_exe=samtools_exe,
            bedgraphtobigwig_exe=bedgraphtobigwig_exe,
            partition_length=partition_length,
            max_readlet_size=max_readlet_size,
            readlet_config_size=readlet_config_size,
            min_readlet_size=min_readlet_size,
            readlet_interval=readlet_interval,
            cap_size_multiplier=cap_size_multiplier,
            max_intron_size=max_intron_size,
            min_intron_size=min_intron_size, min_exon_size=min_exon_size,
            library_size=library_size, 
            search_filter=search_filter,
            motif_search_window_size=motif_search_window_size,
            max_gaps_mismatches=max_gaps_mismatches,
            motif_radius=motif_radius,
            genome_bowtie1_args=genome_bowtie1_args,
            transcriptome_bowtie2_args=transcriptome_bowtie2_args,
            count_multiplier=count_multiplier,
            junction_criteria=junction_criteria,
            indel_criteria=indel_criteria,
            tie_margin=tie_margin,
            normalize_percentile=normalize_percentile,
            transcriptome_indexes_per_sample=transcriptome_indexes_per_sample,
            drop_deletions=drop_deletions,
            do_not_output_bam_by_chr=do_not_output_bam_by_chr,
            do_not_output_ave_bw_by_chr=do_not_output_ave_bw_by_chr,
            do_not_drop_polyA_tails=do_not_drop_polyA_tails,
            deliverables=deliverables, bam_basename=bam_basename,
            bed_basename=bed_basename, tsv_basename=tsv_basename)
        engine_base_checks = {}
        for i in rc.ids:
            engine_base_checks[i] = engine_bases[i].check_program
        if base.check_curl_on_engines:
            apply_async_with_errors(rc, rc.ids, engine_base_checks,
                'curl', 'cURL', '--curl', entered_exe=base.curl_exe,
                reason=base.check_curl_on_engines, is_exe=is_exe, which=which)
        engine_base_checks = {}
        for i in rc.ids:
            engine_base_checks[i] = engine_bases[i].check_s3
        if base.check_s3_on_engines:
            apply_async_with_errors(rc, rc.ids, engine_base_checks,
                reason=base.check_curl_on_engines, is_exe=is_exe, which=which)
        raise_runtime_error(engine_bases)
        print_to_screen(base.detect_message)
        self._json_serial = {}
        step_dir = os.path.join(temp_base_path, 'rna', 'steps')
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
    def __init__(self, manifest, output_dir, isofrag_idx=None,
        intermediate_dir='./intermediate', force=False, aws_exe=None,
        profile='default', region=None,
        service_role=None, instance_profile=None,
        verbose=False, nucleotides_per_input=8000000, gzip_input=True,
        do_not_bin_quals=False, short_read_names=False, skip_bad_records=False,
        bowtie1_exe=None, bowtie_idx='genome', bowtie1_build_exe=None,
        bowtie2_exe=None, bowtie2_build_exe=None, k=1, bowtie2_args='',
        samtools_exe=None, bedgraphtobigwig_exe=None,
        partition_length=5000, max_readlet_size=25,
        readlet_config_size=32, min_readlet_size=15, readlet_interval=4,
        cap_size_multiplier=1.2, max_intron_size=500000, min_intron_size=10,
        min_exon_size=9, library_size=40, search_filter='none',
        motif_search_window_size=1000, max_gaps_mismatches=3,
        motif_radius=5, genome_bowtie1_args='-v 0 -a -m 80',
        transcriptome_bowtie2_args='-k 30', tie_margin=6, count_multiplier=15,
        junction_criteria='0.5,5', indel_criteria='0.5,5',
        normalize_percentile=0.75, transcriptome_indexes_per_sample=500,
        drop_deletions=False, do_not_output_bam_by_chr=False,
        do_not_output_ave_bw_by_chr=False, do_not_drop_polyA_tails=False,
        deliverables='idx,tsv,bed,bw', bam_basename='alignments',
        bed_basename='', tsv_basename='', log_uri=None, ami_version='3.8.0',
        visible_to_all_users=False, tags='', name='Rail-RNA Job Flow',
        action_on_failure='TERMINATE_JOB_FLOW', hadoop_jar=None,
        master_instance_count=1, master_instance_type='c1.xlarge',
        master_instance_bid_price=None, core_instance_count=1,
        core_instance_type=None, core_instance_bid_price=None,
        task_instance_count=0, task_instance_type=None,
        task_instance_bid_price=None, ec2_key_name=None, ec2_subnet_id=None,
        ec2_master_security_group_id=None, ec2_slave_security_group_id=None,
        keep_alive=False, termination_protected=False, check_manifest=True,
        no_direct_copy=False, consistent_view=False,
        intermediate_lifetime=4, max_task_attempts=4, dbgap_key=None,
        secure_stack_name=None):
        base = RailRnaErrors(manifest, output_dir, isofrag_idx=isofrag_idx,
            intermediate_dir=intermediate_dir,
            force=force, aws_exe=aws_exe, profile=profile,
            region=region, service_role=service_role,
            instance_profile=instance_profile, verbose=verbose,
            max_task_attempts=max_task_attempts, dbgap_key=dbgap_key)
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
            ec2_key_name=ec2_key_name, ec2_subnet_id=ec2_subnet_id,
            ec2_master_security_group_id=ec2_master_security_group_id,
            ec2_slave_security_group_id=ec2_slave_security_group_id,
            keep_alive=keep_alive, termination_protected=termination_protected,
            consistent_view=consistent_view,
            no_direct_copy=no_direct_copy,
            intermediate_lifetime=intermediate_lifetime,
            secure_stack_name=secure_stack_name)
        RailRnaPreprocess(base, nucleotides_per_input=nucleotides_per_input,
            gzip_input=gzip_input, do_not_bin_quals=do_not_bin_quals,
            short_read_names=short_read_names,
            skip_bad_records=skip_bad_records)
        RailRnaAlign(base, elastic=True, bowtie1_exe=bowtie1_exe,
            bowtie_idx=bowtie_idx, bowtie1_build_exe=bowtie1_build_exe,
            bowtie2_exe=bowtie2_exe, bowtie2_build_exe=bowtie2_build_exe,
            k=k, bowtie2_args=bowtie2_args, samtools_exe=samtools_exe,
            bedgraphtobigwig_exe=bedgraphtobigwig_exe,
            partition_length=partition_length,
            max_readlet_size=max_readlet_size,
            readlet_config_size=readlet_config_size,
            min_readlet_size=min_readlet_size,
            readlet_interval=readlet_interval,
            cap_size_multiplier=cap_size_multiplier,
            max_intron_size=max_intron_size,
            min_intron_size=min_intron_size,
            search_filter=search_filter,
            min_exon_size=min_exon_size, library_size=library_size, 
            motif_search_window_size=motif_search_window_size,
            max_gaps_mismatches=max_gaps_mismatches,
            motif_radius=motif_radius,
            genome_bowtie1_args=genome_bowtie1_args,
            transcriptome_bowtie2_args=transcriptome_bowtie2_args,
            count_multiplier=count_multiplier,
            junction_criteria=junction_criteria,
            indel_criteria=indel_criteria,
            tie_margin=tie_margin,
            normalize_percentile=normalize_percentile,
            transcriptome_indexes_per_sample=transcriptome_indexes_per_sample,
            drop_deletions=drop_deletions,
            do_not_output_bam_by_chr=do_not_output_bam_by_chr,
            do_not_output_ave_bw_by_chr=do_not_output_ave_bw_by_chr,
            do_not_drop_polyA_tails=do_not_drop_polyA_tails,
            deliverables=deliverables, bam_basename=bam_basename,
            bed_basename=bed_basename, tsv_basename=tsv_basename,
            s3_ansible=ab.S3Ansible(aws_exe=base.aws_exe,
                                        profile=base.profile))
        raise_runtime_error(base)
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
                    base.hadoop_jar, _elastic_step_dir,
                    reducer_count, base.intermediate_dir, unix=True,
                    no_direct_copy=base.no_direct_copy
                ) + \
                steps(
                    RailRnaAlign.protosteps(base, push_dir, elastic=True),
                    base.action_on_failure,
                    base.hadoop_jar, _elastic_step_dir,
                    reducer_count, base.intermediate_dir, unix=True,
                    no_direct_copy=base.no_direct_copy
                )
        self._json_serial['AmiVersion'] = base.ami_version
        self._json_serial['ServiceRole'] = base.service_role
        self._json_serial['JobFlowRole'] = base.instance_profile
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
        self._json_serial['BootstrapActions'] = (
                RailRnaElastic.prebootstrap(base)
                + RailRnaAlign.bootstrap(base)
                + RailRnaPreprocess.srabootstrap(base)
                + RailRnaElastic.bootstrap(base)
            )
        self.base = base
    
    @property
    def json_serial(self):
        return self._json_serial
