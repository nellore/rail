Rail
====

This is the source repository for the Rail-*, which currently includes just
the Rail-RNA suite of RNA-seq analysis software.

Design
------

Here's a flowchart showing how all the pipeline stages relate to each other,
and which source file each is implemented in.  The source files referred to
are in the `src/rnawesome` subdirectory (note: rename this to remove
reference to `rnawesome`).  More detailed descriptions of the
input and output tuples for each step can be found in the source file
headers.

![Design flowchart](doc/design_figure.png "Design flowchart")

Important outputs include:
* `bigBed` (`.bb`) files with per-sample coverage, output by `normalize2.py`
* A `.tsv` file with per-sample normalization factors output by
  `post_normalize.py`
* Some yet-to-be determined format for expressing splicing events, including
  splicing events that co-occurred on one sequencing fragment, output by
  `intron2.py`

Repository layout
-----------------

* `src`: All sources, including Python, R and C
* `data`: Manifest files describing big datasets and some scripts for running
  them through this pipeline as well as Myrna
* `doc`: Documentation
* `emr`: Scripts related to EMR mode
* `example`: Some examples of how to run the pipeline
* `hadoop`: Source related to the modes that use Hadoop
* `lib`: Required Java libraries
* `old`: Deprecated code from pre-Fall-2013
* `package`: Scripts for packaging up the tool and reference archives
* `tools`: Scripts for building some helpful tools

TODO list
---------

* Carefully audit the code in `src/rnawesome/align.py` and add a comprehensive
  suite of unit tests.
* Move this TODO list to GitHub issue tracker

Running the drosophila example
------------------------------

Tool prerequisites:
* [Python]
* [SWIG]
* [Bowtie]
* [bedToBigBed]

[Python]: http://www.python.org
[SWIG]: http://www.swig.org
[Bowtie]: http://bowtie-bio.sourceforge.net/index.shtml
[bedToBigBed]: http://hgdownload.cse.ucsc.edu/admin/exe/

Data prerequisites:
* Expanded dm3 [iGenomes] tarball somewhere on your machine

[iGenomes]: http://support.illumina.com/sequencing/sequencing_software/igenome.ilmn

Configuration:
* Set `TORNADO_HOME` (note: rename!) environment variable to point to base of this checkout
* Set `IGENOMES_HOME` environment variable to point to directory that has
  `Drosophila_melanogaster` from iGenomes tarball as a subdirectory

To prepare for the test:

    cd $TORNADO_HOME
    make -C src
    cd example/dmel_flux
    make
    (alternately, do "make -jN" where N = number of cores you can spare)

That last make command will take a while.  Among other things, it calls Flux
Simulator to simulate some RNA-seq reads.  To run the test do the following
from `$TORNADO_HOME/example/dmel_flux`):

    sh local.sh

Contributors
------------

* [Ben Langmead]
* Jamie Morton
* Abhi Nellore
* [Jeff Leek]
* [Alyssa Frazee]

[Ben Langmead]: http://www.cs.jhu.edu/~langmea/index.shtml
[Jeff Leek]: http://www.biostat.jhsph.edu/~jleek/
[Alyssa Frazee]: http://alyssafrazee.com
