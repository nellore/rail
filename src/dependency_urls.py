"""
dependency_urls.py
Part of Rail-RNA

Defines URLs at which dependencies may be downloaded when installing Rail-RNA.
Increasing version numbers of dependencies should formally increase Rail-RNA's
version.
"""

linux_dependencies = {
    'bowtie1' : ('http://downloads.sourceforge.net/project/bowtie-bio/'
                 'bowtie/1.1.1/bowtie-1.1.1-linux-x86_64.zip'),
    'bowtie2' : ('http://downloads.sourceforge.net/project/bowtie-bio/'
                 'bowtie2/2.2.5/bowtie2-2.2.5-linux-x86_64.zip'),
    'bedgraphtobigwig' : ('http://hgdownload.cse.ucsc.edu/admin/exe/'
                          'linux.x86_64/bedGraphToBigWig'),
    'aws_cli' : 'http://s3.amazonaws.com/aws-cli/awscli-bundle.zip',
    'samtools' : ('http://downloads.sourceforge.net/project/samtools/samtools/'
                  '1.2/samtools-1.2.tar.bz2'),
    'pypy' : ('https://bitbucket.org/squeaky/portable-pypy/downloads/'
              'pypy-2.5-linux_x86_64-portable.tar.bz2'),
}

mac_dependencies = {
    'bowtie1' : ('http://downloads.sourceforge.net/project/bowtie-bio/'
                 'bowtie/1.1.1/bowtie-1.1.1-macos-x86_64.zip'),
    'bowtie2' : ('http://downloads.sourceforge.net/project/bowtie-bio/'
                 'bowtie2/2.2.5/bowtie2-2.2.5-macos-x86_64.zip'),
    'bedgraphtobigwig' : ('http://hgdownload.cse.ucsc.edu/admin/exe/'
                          'macOSX.x86_64/bedGraphToBigWig'),
    'aws_cli' : 'http://s3.amazonaws.com/aws-cli/awscli-bundle.zip',
    'samtools' : ('http://downloads.sourceforge.net/project/samtools/samtools/'
                  '1.2/samtools-1.2.tar.bz2'),
    'pypy' : ('http://bitbucket.org/pypy/pypy/downloads/'
              'pypy-2.5.0-osx64.tar.bz2')
}