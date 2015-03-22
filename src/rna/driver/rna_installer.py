#!/usr/bin/env python
"""
rna_installer.py
Part of Rail-RNA

Contains a class for installing Rail-RNA.
"""
import sys
import contextlib
import os
from rna_config import *
base_path = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
import site
site.addsitedir(base_path)
import dependency_urls
from distutils.util import strtobool
from dooplicity.tools import which, register_cleanup
import zipfile
import shutil
import subprocess
from version import version_number
import multiprocessing
import tempfile

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
        self.zip_name = os.path.abspath(zip_name)
        self.curl_exe = curl_exe
        log_dir = tempfile.mkdtemp()
        self.log_file = os.path.join(log_dir, 'rail-rna_install.log')
        self.log_stream = open(log_file, 'w')
        register_cleanup(log_dir)

    def __enter__(self):
        return self

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
        print_to_screen('Downloading %s...' % name)
        command = [self.curl_exe, '-L', '-O', url]
        filename = url.partition('/')[2]
        try:
            subprocess.check_output(command, stderr=sys.stderr)
        except subprocess.CalledProcessError as e:
            print_to_screen(
                    ('Error encountered downloading file %s; exit '
                     'code was %d; command invoked was "%s".') %
                        (url, e.returncode, ' '.join(command))
                )
            print_to_screen('Make sure web access is available.')
            self._bail()
        else:
            # Explode
            explode_command = None
            if url[:-8] == '.tar.bz2':
                explode_command = ['tar', 'xvjf', filename]
            elif url[:-7] == '.tar.gz' or url[:-4] == '.tgz':
                explode_command = ['tar', 'xvjf', filename]
            elif url[:-4] == '.zip':
                try:
                    zipfile.extractall(filename)
                except Exception as e:
                    print_to_screen('Error encountered exploding %s.'
                                        % filename)
                    self._bail()
                finally:
                    os.remove(filename)
            if explode_command is not None:
                try:
                    subprocess.check_output(explode_command,
                                            stderr=self.log_stream)
                except subprocess.CalledProcessError as e:
                    print_to_screen(
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
                self._bail()
        if self._yes_no_query(
                'Rail-RNA can be installed for all users or for just the '
                'current user. Install for all users?'
            ):
            if os.getuid():
                print_to_screen('Rerun with sudo privileges to install '
                                'for all users.')
                sys.exit(0)
            install_dir = '/usr/local'
            exe_dir = '/usr/local/bin'
        else:
            install_dir = os.path.expanduser('~/.local/bin')
            exe_dir = None
        full_install_dir = os.path.join(install_dir, 'rail-rna')
        if os.path.exists(full_install_dir):
            if self._yes_no_query(
                    'An installation exists at %s. Overwrite?'
                    % full_install_dir
                ):
                try:
                    shutil.rmtree(full_install_dir)
                except OSError:
                    # Handle this if directory creation fails
                    pass
            else:
                print_to_screen('Exiting.')
                sys.exit(0)
        print_to_screen('Extracting Rail-RNA...')
        try:
            os.makedirs(full_install_dir)
        except OSError:
            print_to_screen(('Problem encountered trying to create '
                             'directory %s for installation.'))
            self._bail()
        with cd(full_install_dir):
            zipfile.extractall(zip_name)
            print_to_screen('Installing dependencies.')
            self._grab_and_explode(self.depends['bowtie1'], 'Bowtie 1')
            self._grab_and_explode(self.depends['bowtie2'], 'Bowtie 2')
            self._grab_and_explode(self.depends['bedgraphtobigwig'],
                                    'BedGraphToBigWig')
            self._grab_and_explode(self.depends['pypy'], 'PyPy')
            self._grab_and_explode(self.depends['samtools'], 'SAMTools')
        # Have to make SAMTools (annoying; maybe change this)
        samtools_dir = os.path.join(full_install_dir,
                self.depends['samtools'].partition('/')[2].split('.')[0]
            )
        with cd(samtools_dir):
            # Make on all but one cylinder
            thread_count = max(1, multiprocessing.cpu_count() - 1)
            samtools_command = ['make', '-j%d' % thread_count]
            try:
                subprocess.check_output(['make', '-j%d' % thread_count],
                                            stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e
                print_to_screen(
                        ('Error encountered making SAMTools; exit '
                         'code was %d; command invoked was "%s".') %
                            (e.returncode, ' '.join(samtools_command))
                    )
                self._bail()

    def uninstall(self):
        """ Uninstalls Rail-RNA. """

    def __exit__(self, type, value, traceback):
        try:
            self.log_stream.close()
        except:
            pass

