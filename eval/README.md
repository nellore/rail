Note on preprint versioning
-----
This version of the readme is describes how to reproduce results in v2 of the [Rail-RNA preprint](http://biorxiv.org/content/early/2015/05/07/019067). To view the README that reproduces the results contained in v1, see [this](https://github.com/nellore/rail/tree/7db1e4e506224e156ae1ac271ceb88d5e465b7b9/eval) commit.

Reproducing preprint results
-----
**Use Rail-RNA v0.1.9b**. Rail-RNA has several dependencies. In the experiments we conducted, Rail-RNA wrapped Bowtie 1 v1.1.1, Bowtie 2 v2.2.5, and SAMTools v1.2, and PyPy 2.5. Version 2.1.0 of TopHat 2 was used and, like Rail-RNA, it wrapped Bowtie 2 v2.2.5. Version 2.4.2a of STAR was used, and version 0.1.6-beta of HISAT was used. Subjunc is a
tool in the Subread package, and version 1.4.6-p4 of the Subread package was used. Flux Simulator 1.2.1 was used to obtain simulated samples as described below.

To reproduce results from v2 of the Rail-RNA [preprint](http://biorxiv.org/content/early/2015/05/07/019067), perform the following steps. Note that input and output directories in scripts may need to be changed. If you're having trouble, ask questions in our [Gitter](https://gitter.im/nellore/rail).

Generating genome indexes for all experiments
-----
Use `create_indexes.sh`. Refer to comments in that script for further instructions.

Scaling experiments
-----
1. Start with the Myrna-style manifest file GEUVADIS_all_samples.manifest composed of the URLs of all GEUVADIS samples (minus one; see note under the section All-of-GEUVADIS run below). Invoke
    ```
    python geuvadis_lab_ethnicity_sex_all.py
    ```
to create from it another manifest file GEUVADIS_all_descriptive.manifest with more descriptive sample names including gender, ethnicity, and lab. The script uses data made public by the GEUVADIS consortium and included in this directory. Details on data sources are in the script.
2. Run

    ```
    python geuvadis_lab_ethnicity_sex.py
    ```
to generate a random sample A of the GEUVADIS corpus with 112 paired-end datasets, 16 from each lab, a random sample B of A
with 56 datasets, 8 from each lab, and a random sample C of B with 28 datasets, 4 from each lab. (Again, see script for data sources.) The corresponding manifest files are stored in `GEUVADIS_112.manifest`, `GEUVADIS_56.manifest`, and `GEUVADIS_28.manifest`. Optionally, use `make_sample_labels_consistent.py` to ensure that sample labels from all manifest files are the same for the same datasets. (Otherwise, replicate numbers could differ between two manifests for the same dataset.)
3. Change the output S3 bucket and home directory of Rail-RNA's source code (`rail/src`) in `prep_scaling_experiments.sh`. Run

    ```
    sh prep_scaling_experiments.sh
    ```
to submit all jobs that preprocess GEUVADIS data for scaling experiments to Elastic MapReduce (EMR). This requires an account with [Amazon Web Services](http://aws.amazon.com/), whose charges will apply. Note also that `--ec2-key-name` in all commands from `prep_scaling_experiments.sh` should either be removed or changed to a valid EC2 key name. The key name permits SSHing to EC2 instances while the preprocessing jobs are running. See (this)[http://docs.aws.amazon.com/gettingstarted/latest/wah/getting-started-prereq.html] page for more information.

4. Again, change the output S3 bucket and home directory of Rail-RNA's source code in `submit_scaling_experiments.sh`, and run

    ```
    sh submit_scaling_experiments.sh
    ```
to execute scaling experiments on EMR. The 28 GEUVADIS datasets C are analyzed on three (virtualized) computer clusters of different sizes: 10, 20, and 40 c3.2xlarge instances. The 56 datasets B and 128 datasets A are also analyzed on 40 c3.2xlarge instances. The experiments provide information on how well Rail-RNA scales with respect to both input data volume and cluster size. For the scalability runs discussed in the preprint, wall clock times were read off the EMR interface. See `scaling_experiment_results.tsv` for the results. They are depicted in Figure 1 of the preprint. To regenerate the plots in Figure 1, run `performance.nb` in [Wolfram Mathematica](http://www.wolfram.com/mathematica/) >= 9.

Accuracy experiments
-----
Prerequisite: manifest file `GEUVADIS_112.manifest` from the scaling experiments above.

1. Run

    ```
    python generate_bioreps.py --help
    ```
to view instructions on how to regenerate the 112 simulated bioreplicates whose coverage distributions mimic the samples from `GEUVADIS_112.manifest`. [Flux Simulator](http://sammeth.net/confluence/display/SIM/Home) 1.2.1 is required. The script uses Flux Simulator's built-in 76-bp error model to generate reads with substitution errors. See the docstring in `get_error_distribution.py` for information on how to reproduce the statements on substitution error rates in the supplement of v2 of the [Rail-RNA preprint](http://biorxiv.org/content/early/2015/05/07/019067).
2. In v2 of the [Rail-RNA preprint]([Rail-RNA preprint](http://biorxiv.org/content/early/2015/05/07/019067)), 20 samples are selected at random the 112 simulated samples, and various alignment protocols are executed on them. These protocols use [TopHat 2](http://ccb.jhu.edu/software/tophat/index.shtml), [STAR](https://github.com/alexdobin/STAR), [HISAT](http://ccb.jhu.edu/software/hisat/index.shtml), and [Subjunc](http://subread.sourceforge.net/). Use the script `create_single_sample_sim_commands.py` to generate a set of commands that output the results of these protocols to the same directory, which can be toggled with the `--output-dir` command-line parameter. `create_single_sample_sim_commands.py` depends on five scripts also contained in this directory: `run_single_sample_{star,tophat,hisat,rail,subjunc}_sim.sh`. These five scripts invoke three other scripts to compute accuracy metrocs: `spliced_read_recovery_performance.py`, `intron_recovery_performance.py`, and `mapping_accuracy.py`. See the docstrings of these Python scripts for information on what precisely they measure. For an example of how to use `create_single_sample_sim_commands.py`, see `create_commands_for_all_sample_sims.sh`, which we used locally to generate some commands for our experiments.
3. Rail-RNA borrows strength from across samples by accumulating a master list of introns for realigning reads. Use `run_112_sim.sh` to run Rail-RNA on the 112 simulated bioreps on EMR in two ways: (a) with an intron filter that, before realignment, eliminates introns occurring in no more than 5% of samples and initially detected in no more than 5 reads per sample; (b) without an intron filter. See the comments in that script for information on command-line parameters; note that it stages the FASTQs generated by `generate_bioreps.py` on S3.
4. After the experiments are complete, use `grab_112_sim_results.sh` to download results for the same twenty simulated samples selected at random in `create_single_sample_sim_commands.py` and compute accuracy metrics. These results are influenced by the presence of the other samples, unlike the results of the commands generated by `create_single_sample_sim_commands.py`. Check the comments of `grab_112_sim_results.sh` to learn about command-line parameters; in particular, the output directory of the experiments in Step 3 and the path to the 112 simulated bioreplicates should be specified. The script uses `count_introns.py` to get information on canonical and noncanonical splice sites on the filtered and unfiltered intron lists. Figure 2b from the preprint depicts intron call precision/recall and canonical/noncanonical status written to `withoutfilter/count_introns` (before intron filtering) and `withfilter/count_introns` (after intron filtering) by `grab_112_sim_results.sh`.
5. Use `harvest_accuracies.py` to construct tables of precisions, recalls, and F-scores as well as their means and standard deviations across the twenty samples for the various accuracy metrics used in the preprint. The command-line parameters `--sam-dir` and `--aligners` must be specified; the former is the directory where either `grab_112_sim_results.sh` or the commands generated by `create_single_sample_sim_commands.py` dumped their output, while the latter should be a space-separated list of aligners from `{star,tophat,hisat,rail,subjunc}` whose output is contained in the `--sam-dir`. In particular, `--aligner` must be only `rail` when `--sam-dir` is the output directory of `grab_112_sim_results.sh`, and `--aligner` should be `star tophat hisat rail subjunc` when `--sam-dir` is the output directory from `create_single_sample_sims.py`. The tables output by `harvest_accuracies.py` are reproduced in the supplement of the Rail-RNA preprint. Means and standard deviations are depicted in Figure 2a.
6. Apply `consolidate_performances.sh` to copy all performance results of `grab_112_sim_results.sh` and the commands generated by `create_single_sample_sims.py` to the same directory. See the comments in `consolidate_performances.sh` for help making it from this step to the next.
7. Open `performance.nb` in [Wolfram Mathematica](http://www.wolfram.com/mathematica/) >=9 to read this performance data and generate Figures 7 and 8 from the preprint. Note that wherever plots are unannotated, Keynote was used to add back annotations.
8. Use `harvest_coverages_from_performance.sh` to obtain tables that give spliced alignment accuracy metrics at various read coverages of introns. This data is used to generate Figures 5 and 6 from the preprint, recording mean values and standard deviations. The full tables are provided in the supplement.

All-of-GEUVADIS run
-----
The file GEUVADIS_all_samples.manifest from which 112 samples are obtained for scaling experiments is missing one GEUVADIS sample because its FASTQ could not be found on the server. We ultimately realized that the missing sample was paired-end and divided into two FASTQs on the server. The original manifest file with the missing sample `GEUVADIS_all_samples.manifest` has been preserved so manifest files for scaling and accuracy experiments may be regenerated, but we use a different manifest file for our all-of-GEUVADIS run. To generate it:

1. Start with the Myrna-style manifest file GEUVADIS_all_samples_revised.manifest composed of the URLs of all GEUVADIS samples. The missing sample from GEUVADIS_all_samples.manifest is added to the top. Invoke
    ```
    python geuvadis_lab_ethnicity_sex_all_revised.py
    ```
2. Change the output bucket and argument of `--ec2-key-name` in `preprocess_all_of_GEUVADIS.sh`, and run the preprocess job flow contained in that script.
3. Change the output bucket and argument of `--ec2-key-name` in `submit_all_of_geuvadis_job.sh`, and run the alignment job flow contained in that script.

`--ec2-key-name` may also be removed in the scripts above. This is a PEM file the user generates to be able to SSH to an Elastic MapReduce cluster. See (this)[http://docs.aws.amazon.com/gettingstarted/latest/wah/getting-started-prereq.html] page for more information.

Junction comparison analysis
----
1. Download all junction BEDs form the all-of-GEUVADIS run and place them in the same directory `D`.
2. Run `sh generate_geuvadis_intron_table.sh D`, where `D` is the directory from step 1. Two files are output: one (`intron_table.tsv`) is a list of introns overlapped by Rail-RNA's alignments across all of GEUVADIS; the other (`GEUVADIS_introns_v4.04.2015.tsv.gz`) is a matrix whose (i, j)th element is the number of reads overlapping intron i in sample j.
3. In `junction_comparison.R`, under `## read in Rail jxns`, write the correct paths to `intron_table.tsv` and `GEUVADIS_introns_v4.04.2015.tsv.gz`.
4. Run `junction_comparison.R`.

derfinder analysis
----

1. Use `railDER/run-all.sh` to load the coverage from the bigWig files and identify the Expressed Regions (ERs). To run it, use:

    ```
    cd railDER
    sh run-all.sh railGEU
    ```
This will generate the files in `railDER/railGEU/CoverageInfo` and `railDER/railGEU/regionMatrix`. `railDER/run-all.sh` is based on `derSoftware/run-all.sh` described at [derSoftware](http://leekgroup.github.io/derSoftware/) which is the supplementary website for the `derfinder` paper.
2. Run `railGEU/fixSampleNames/fixSampleNames.R` to match the sample names from the GEUVADIS `ballgown` object and the names used by Rail-RNA in the manifest files.

    ```
    cd railDER/railGEU
    Rscript fixSampleNames.R
    ```
This will generate the `GRanges` object with the ERs in `railDER/railGEU/fixSampleNames/regions.Rdata` and the coverage matrix `railDER/railGEU/fixSampleNames/coverageMatrix.Rdata`. The coverage matrix file is available via Figshare at [dx.doi.org/10.6084/m9.figshare.1405638](http://dx.doi.org/10.6084/m9.figshare.1405638).
3. `derfinder_analysis.R` analyzes the resulting ERs and generates several plots included in the preprint and supplementary material. Use

    ```
    Rscript derfinder_analysis.R
    ```
Alternatively, submit a SGE cluster job using:

    ```
    qsub derfinder_analysis.sh
    ```
This generates `ssMat_geuvadis.rda` that has for all the 290416 ERs the percent of variance explained by population, each of the 14 technical variables and the residual variation: population, RIN value, RNA extraction batch, RNA concentration, RNA quantity used, Library preparation date, primer index, method concentration measure, library concentration, library size, library concentration used, cluster kit, sequencing kit, cluster density, lane, and residual variation. It then creates boxplots saved in `r2_boxplots_overall.pdf` and `r2_boxplots_overall_noPop.pdf` (without population variable) as well as a venn diagram of the ERs overlapping annotation features from ENSEMBL v75 saved in `ensemblVenn.pdf`. Overlap information is saved in `ensemblAnno.Rdata` while explained variance is saved in `ssMat_geuvadis.rda` and `ssMat_geuvadis_noPop.rda`. `derfinder_analysis.o4348498` and `derfinder_analysis.e4348498` are the corresponding log files.

GEUVADIS read count histogram (Figure 5 from preprint)
----
Download the 667 paired-end GEUVADIS sample FASTQ.gzs contained in `GEUVADIS_all_descriptive_revised.manifest`. Run the command `gzip -cd | wc -l` on each file ending with `_1.fastq.gz` from the manifest, and divide the result by 2 to obtain the number of reads in the corresponding sample. We plotted the histogram in [Wolfram Mathematica](http://www.wolfram.com/mathematica/) v10 using `performance.nb`.

Exon-exon junction count histogram (Figure 6 from preprint)
----
Prerequisite: Step 1 under Accuracy experiments above.

See the docstring of `intron_histogram.py` for information on how to regenerate the histogram data. This requires that We plotted the histogram in [Wolfram Mathematica](http://www.wolfram.com/mathematica/) v10 using `performance.nb`.
