#!/usr/bin/env python
"""
rna_installer.py
Part of Rail-RNA

Contains a class for installing Rail-RNA.
"""
import sys
import contextlib
import os
base_path = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
utils_path = os.path.join(base_path, 'rna', 'utils')
import site
site.addsitedir(base_path)
site.addsitedir(utils_path)
import dependency_urls
from distutils.util import strtobool
from dooplicity.tools import which, register_cleanup
import zipfile
import shutil
import subprocess
from version import version_number
import multiprocessing
import tempfile
from tempdel import remove_temporary_directories

@contextlib.contextmanager
def cd(dir_name):
    """ Changes directory in a context only. Borrowed from AWS CLI code. """
    original_dir = os.getcwd()
    os.chdir(dir_name)
    try:
        yield
    finally:
        os.chdir(original_dir)

class RailRnaInstaller(object):
    """ Installs Rail-RNA and its assorted dependencies.

        Init vars
        -------------
        archive_name: path to (currently executing) zip containing Rail-RNA
    """

    def __init__(self, zip_name, curl_exe=None):
        print_to_screen('{0} Rail-RNA v{1} Installer'.format(
                                        , version_number)
                                    )
        if sys.platform in ['linux', 'linux2']:
            self.depends = dependency_urls.linux_dependencies
            self.linux = True
        elif sys.platform == 'darwin':
            self.depends = dependency_urls.mac_dependencies
            self.linux = False
        else:
            print_to_screen(
                    'Rail-RNA cannot be installed because it is not supported '
                    'by your OS. Currently supported OSes are Mac OS X and '
                    'Linux.'
                )
            sys.exit(1)
        self.zip_name = os.path.abspath(zip_name)
        self.curl_exe = curl_exe
        log_dir = tempfile.mkdtemp()
        self.log_file = os.path.join(log_dir, 'rail-rna_install.log')
        self.log_stream = open(log_file, 'w')
        register_cleanup(remove_temporary_directories, [log_dir])

    def __enter__(self):
        return self

    def _print_to_screen_and_log(message, **kwargs):
        print >>self.log_stream, message
        print_to_screen(message, **kwargs)

    def _bail(self):
        """ Copy log to some temporary dir and GTFO. """
        new_log_file = os.path.join(tempfile.mkdtemp(), 'rail-rna_install_log')
        shutil.copyfile(self.log_file, new_log_file)
        print_to_screen('Installation log may be found at %s.' % new_log_file)
        sys.exit(1)

    def _yes_no_query(self, question):
        """ Gets a yes/no answer from the user.

            question: string with question to be printed to console

            Return value: boolean
        """
        while True:
            sys.stdout.write('%s [y/n]: ' % question)
            try:
                return strtobool(raw_input().lower())
            except ValueError:
                sys.stdout.write('\nPlease enter \'y\' or \'n\'.\n')

    def _grab_and_explode(self, url, name):
        """ Special method for grabbing and exploding a package, if necessary.

            Does not verify URL since these are preverified. Package is
            downloaded to current directory.

            url: url to grab
            name: name of download

            No return value
        """
        self._print_to_screen_and_log('Downloading %s...' % name)
        command = [self.curl_exe, '-L', '-O', url]
        filename = url.partition('/')[2]
        try:
            subprocess.check_output(command, stderr=self.log_stream)
        except subprocess.CalledProcessError as e:
            self._print_to_screen_and_log(
                    ('Error encountered downloading file %s; exit '
                     'code was %d; command invoked was "%s".') %
                        (url, e.returncode, ' '.join(command))
                )
            self._print_to_screen_and_log('Make sure web access is available.')
            self._bail()
        else:
            # Explode
            self._print_to_screen_and_log('Exploding %s...' % name)
            explode_command = None
            if url[:-8] == '.tar.bz2':
                explode_command = ['tar', 'xvjf', filename]
            elif url[:-7] == '.tar.gz' or url[:-4] == '.tgz':
                explode_command = ['tar', 'xvjf', filename]
            elif url[:-4] == '.zip':
                try:
                    zipfile.extractall(filename)
                except Exception as e:
                    self._print_to_screen_and_log(
                            'Error encountered exploding %s.' % filename
                        )
                    self._bail()
                finally:
                    os.remove(filename)
            if explode_command is not None:
                try:
                    subprocess.check_output(explode_command,
                                            stderr=self.log_stream)
                except subprocess.CalledProcessError as e:
                    self._print_to_screen_and_log(
                        ('Error encountered exploding file %s; exit '
                         'code was %d; command invoked was "%s".') %
                            (filename, e.returncode, ' '.join(explode_command))
                    )
                    self._bail()
                finally:
                    os.remove(filename)

    def install(self):
        """ Installs Rail-RNA and all its dependencies. """
        if self.curl_exe is None:
            self.curl_exe = which('curl')
            if self.curl_exe is None:
                print_to_screen('Rail-RNA\'s installer requires Curl. '
                                'Download it at '
                                'http://curl.haxx.se/download.html.')
                sys.exit(1)
        if self._yes_no_query(
                'Rail-RNA can be installed for all users or for just the '
                'current user. Install for all users?'
            ):
            if os.getuid():
                print_to_screen('Rerun with sudo privileges to install '
                                'for all users.')
                sys.exit(0)
            install_dir = '/usr/local'
            local = False
        else:
            install_dir = os.path.abspath(os.path.expanduser('~/.local'))
            local = True
        rail_exe = os.path.join(install_dir, 'bin', 'rail-rna')
        final_install_dir = os.path.join(install_dir, 'rail-rna')
        # Install in a temporary directory first, then move to final dest
        temp_install_dir = tempfile.mkdtemp()
        register_cleanup(remove_temporary_directories, [temp_install_dir])
        if os.path.exists(final_install_dir):
            if self._yes_no_query(
                    'An installation exists at %s. Overwrite?'
                    % final_install_dir
                ):
                try:
                    shutil.rmtree(final_install_dir)
                except OSError:
                    # Handle this if directory creation fails
                    pass
            else:
                print_to_screen('Exiting.')
                sys.exit(0)
        self._print_to_screen_and_log('Extracting Rail-RNA...')
        try:
            os.makedirs(final_install_dir)
        except OSError:
            self._print_to_screen_and_log(
                            ('Problem encountered trying to create '
                             'directory %s for installation.')
                        )
            self._bail()
        else:
            # So it's possible to move temp installation dir there
            os.rmdir(final_install_dir)
            pass
        with cd(temp_install_dir):
            zipfile.extractall(zip_name)
            self._print_to_screen_and_log('Installing dependencies.')
            self._grab_and_explode(self.depends['bowtie1'], 'Bowtie 1')
            self._grab_and_explode(self.depends['bowtie2'], 'Bowtie 2')
            self._grab_and_explode(self.depends['bedgraphtobigwig'],
                                    'BedGraphToBigWig')
            self._grab_and_explode(self.depends['pypy'], 'PyPy')
            self._grab_and_explode(self.depends['samtools'], 'SAMTools')
        # Have to make SAMTools (annoying; maybe change this)
        samtools_dir = os.path.join(temp_install_dir,
                self.depends['samtools'].partition('/')[2].split('.')[0]
            )
        with cd(samtools_dir):
            # Make on all but one cylinder
            thread_count = max(1, multiprocessing.cpu_count() - 1)
            samtools_command = ['make', '-j%d' % thread_count]
            self._print_to_screen_and_log('Making SAMTools...')
            try:
                subprocess.check_output(samtools_command,
                                            stderr=self.log_stream)
            except subprocess.CalledProcessError as e
                self._print_to_screen_and_log(
                        ('Error encountered making SAMTools; exit '
                         'code was %d; command invoked was "%s".') %
                            (e.returncode, ' '.join(samtools_command))
                    )
                self._bail()
        samtools = os.path.join(samtools_dir, 'samtools')
        bowtie1 = os.path.join(temp_install_dir,
                                self.depends['bowtie1'][:-8],
                                'bowtie')
        bowtie1_build = os.path.join(temp_install_dir,
                                self.depends['bowtie1'][:-8],
                                'bowtie-build')
        bowtie2 = os.path.join(temp_install_dir,
                                self.depends['bowtie2'][:-8],
                                'bowtie2')
        bowtie2_build = os.path.join(temp_install_dir,
                                self.depends['bowtie2'][:-8],
                                'bowtie2-build')
        pypy = os.path.join(temp_install_dir,
                                self.depends['pypy'][:-8], 'bin', 'pypy'
                        )
        bedgraphtobigwig = os.path.join(bedgraphtobigwig, 'bedGraphToBigWig')
        os.renames(os.path.join(temp_install_dir, 'bedGraphToBigWig'),
                        bedgraphtobigwig)
        # Write paths to exe_paths
        with open(
                        os.path.join(temp_install_dir, 'exe_paths.py'), 'w'
                    ) as exe_paths_stream:
            print >>exe_paths_stream, (
"""
\"""
exe_paths.py
Part of Rail-RNA

Defines default paths of Rail-RNA's executable dependencies. Set a given
variable equal to None if the default path should be in PATH.
\"""

pypy = {pypy}
aws = None
curl = None
sort = None
bowtie1 = {bowtie1}
bowtie1_build = {bowtie1_build}
bowtie2 = {bowtie2}
bowtie2_build = {bowtie2_build}
samtools = {samtools}
bedgraphtobigwig = {bedgraphtobigwig}
"""
            ).format(pypy=pypy, bowtie1=bowtie1, bowtie1_build=bowtie1_build,
                        bowtie2=bowtie2, bowtie2_build=bowtie2_build,
                        samtools=samtools, bedgraphtobigwig=bedgraphtobigwig)
        # Move to final directory
        try:
            os.renames(temp_install_dir, final_install_dir)
        except OSError:
            self._print_to_screen_and_log(('Problem encountered moving '
                                           'temporary installation %s to '
                                           'final destination %s.') % (
                                                temp_install_dir,
                                                final_install_dir
                                            ))
            self._bail()
        # Create shell-script executable
        with open(rail_exe, 'w') as rail_exe_stream:
            print >>rail_exe_stream, (
"""
#!/usr/bin/env bash

{python_executable} {install_dir} \$@
"""
                ).format(python_executable=sys.executable,
                            install_dir=final_install_dir)
        if local:
            # Have to add Rail to PATH
        self._print_to_screen_and_log(
            'Rail-RNA installation complete.')
        if not which('aws') and self._yes_to_query(
                'AWS CLI is not installed but required for Rail-RNA to work '
                'in its "elastic" mode, on Amazon Elastic MapReduce. '
                'Install now?'
            )
            temp_aws_install_dir = tempfile.mkdtemp()
            register_cleanup(remove_temporary_directories,
                                [temp_aws_install_dir])
            self._print_to_screen_and_log('Installing AWS CLI...')
            with cd(temp_aws_install_dir):
                self._grab_and_explode(self.depends['aws'], 'AWS CLI')
                if local:
                    # Local install
                    aws_command = ['./awscli-bundle/install', '-b',
                                    os.path.abspath(
                                            os.path.expanduser('~/bin/aws')
                                        )]
                else:
                    # All users
                    aws_command = ['./awscli-bundle/install', '-i',
                                '/usr/local/aws', '-b', '/usr/local/bin/aws']
                try:
                    subprocess.check_output(aws_command,
                                                stderr=self.log_stream)
                except subprocess.CalledProcessError as e
                    self._print_to_screen_and_log(
                            ('Error encountered installing AWS CLI; exit '
                             'code was %d; command invoked was "%s".') %
                                (e.returncode, ' '.join(aws_command))
                        )
                    self._bail()
            print_to_screen('Configure the AWS CLI by running '
                            '"aws configure".')
        else:
            print_to_screen('Visit http://docs.aws.amazon.com/cli/latest/'
                            'userguide/installing.html to install the '
                            'AWS CLI later.')
        print_to_screen('Start using Rail by entering "rail-rna".')

    def __exit__(self, type, value, traceback):
        try:
            self.log_stream.close()
        except:
            pass
