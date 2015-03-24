Rail
====
This is the official repo for Rail-*, which currently includes just the Rail-RNA pipeline for spliced alignment of RNA-seq data. Check `releases/` for releases.

Get interested
-----
Rail-RNA's distinguishing features are
* Scalability: built on MapReduce, the software scales to analyze hundreds of RNA-seq samples at the same time. Moreover, through reduced-redundancy analysis, the end-to-end analysis time per sample *decreases* as the number of samples increases.
* Integrative analysis: the software borrows strength across replicates to achieve more accurate splice junction detection, especially at low levels of coverage.
* Mode agnosticism: Rail-RNA integrates its own parallel abstraction layer that allows it to be run in various distributed computing environments, including the Amazon Web Services (AWS) [Elastic MapReduce (EMR) service](http://aws.amazon.com/elasticmapreduce/), or any distributed environment supported by [IPython](http://ipython.org/), including clusters using batch schedulers like PBS or SGE, Message Passing Interface (MPI), or any cluster with a shared filesystem and mutual SSH access. Alternately, Rail-RNA can be run on a single multi-core computer, without the aid of a batch system or MapReduce implementation.

Outputs include
* Alignment BAMs
* Genome coverage bigWigs
* [TopHat](http://ccb.jhu.edu/software/tophat/index.shtml)-like indel and splice junction BEDs.

Read the [preprint](http://www.shouldreallyfinishthis.com/) for more details.

Install
-----
Start with a recent (>= 2009) Mac OS or Linux box. Download [`installers/rail-rna-1.0.0_installer`](https://github.com/buci/rail/blob/master/releases/install_rail-rna-1.0.0?raw=true), change to the directory containing it, and run
```
sudo ./install_rail-rna-1.0.0
```
to install for all users or
```
./install_rail-rna-1.0.0
```
to install for just you. Check
```
./install_rail-rna-1.0.0 -h
```
for more installation options. If the executable doesn't work, you may need [Python](http://www.python.org).

Use
-----
For now, enter
```
rail-rna
```
and follow the instructions. A manual will be posted shortly.

Contributors
-----
* [Abhi Nellore]
* Mark Pritt
* Leo Collado-Torres
* Jamie Morton
* [Jeff Leek]
* [Ben Langmead]

[Abhi Nellore]: https://scholar.google.com/citations?user=XxPWj5oAAAAJ&hl=en
[Leo Collado-Torres]: https://github.com/lcolladotor
[Ben Langmead]: http://www.cs.jhu.edu/~langmea/index.shtml
[Jeff Leek]: http://www.biostat.jhsph.edu/~jleek/
