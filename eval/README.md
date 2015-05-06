Reproducing preprint results
-----
**Use Rail-RNA v0.1.0a.** Rail-RNA has several dependencies. In the experiments we conducted, Rail-RNA wrapped Bowtie 1 v1.1.1, Bowtie 2 v2.2.4, and SAMTools v0.1.19. When run in its `elastic` mode as described below, Rail-RNA used PyPy v2.2.1. In all other cases, Rail-RNA used PyPy v2.4.0. Version 2.0.12 of TopHat 2 was used and, like Rail-RNA, it wrapped Bowtie 2 v2.2.4 and SAMTools v0.1.19. Version 2.4.0j of STAR was used, and version 0.1.5-beta of HISAT was used. Flux Simulator 1.2.1 was used to obtain simulated samples as described below.

To reproduce results from the [preprint](http://finishit.com), perform the following steps. Note that input and output directories in scripts may need to be changed.

Generating genome indexes for all experiments
-----
Use `create_indexes.sh`. Refer to comments in that script for further instructions.

Scaling experiments
-----
1. Start with the Myrna-style manifest file GEUVADIS_all_samples.manifest composed of the URLs of all GEUVADIS samples. Invoke
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
to submit all jobs that preprocess GEUVADIS data for scaling experiments to Elastic MapReduce (EMR). This requires an account with [Amazon Web Services](http://aws.amazon.com/), whose charges will apply. Note also that `--ec2-key-name` in all commands from `prep_scaling_experiments.sh` should either be removed to changed to a valid EC2 key name. The key name permits SSHing to EC2 instances while the preprocessing jobs are running.
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
to view instructions on how to regenerate the 112 simulated bioreplicates whose coverage distributions mimic the samples from `GEUVADIS_112.manifest`. [Flux Simulator](http://sammeth.net/confluence/display/SIM/Home) 1.2.1 is required. 
2. Use the script `run_all_small_data_sims_locally.sh` to run [TopHat 2](http://ccb.jhu.edu/software/tophat/index.shtml), [STAR](https://github.com/alexdobin/STAR), [HISAT](http://ccb.jhu.edu/software/hisat/index.shtml), and Rail-RNA on two samples from the 112 simulated bioreplicates independently. Comments in that script specify precisely its command-line parameters and the versions of all required software used for the preprint. The script invokes two other scripts to compute accuracy metrics: `spliced_read_recovery_performance.py` and `intron_recovery_performance.py`. See the docstrings of these Python scripts for information on what precisely they measure.
3. Rail-RNA borrows strength from across samples by accumulating a master list of introns for realigning reads. Use `run_112_sim.sh` to run Rail-RNA on the 112 simulated bioreps on EMR in two ways: (a) with an intron filter that, before realignment, eliminates introns occurring in no more than 5% of samples and initially detected in no more than 5 reads per sample; (b) without an intron filter. See the comments in that script for information on command-line parameters; note that it stages the FASTQs generated by `generate_bioreps.py` on S3.
4. After the experiments are complete, use `grab_112_sim_results.sh` to download results for the same two samples analyzed in `run_all_small_data_sims_locally.sh` and compute accuracy. These results are influenced by the presence of the other samples, unlike the results of `run_all_small_data_sims_locally.sh`. Check the comments of `grab_112_sim_results.sh` to learn about command-line parameters; in particular, the output directory of the experiments in Step 3 and the path to the 112 simulated bioreplicates should be specified. The script uses `count_introns.py` to get information on canonical and noncanonical splice sites on the filtered and unfiltered intron lists. Figure 2b from the preprint depicts intron call precision/recall and canonical/noncanonical status written to `withoutfilter/count_introns` (before intron filtering) and `withfilter/count_introns` (after intron filtering) by `grab_112_sim_results.sh`.
5. Apply `consolidate_performances.sh` to copy all performance results from `run_all_small_data_sims_locally.sh` and `grab_112_sim_results.sh` in the same directory. Filenames are assigned according to which results correspond to which experiments. The files with `perform_summary` in their filenames contain the precisions and recalls depicted in Figure 2a of the preprint. F-scores were computed with a calculator from these precisions and recalls.
6. Use `harvest_coverages_from_performance.sh` to obtain tables that give spliced alignment accuracy metrics at various read coverages of introns. This data is used to generate Figures 5 and 6 from the preprint.
7. Open `performance.nb` in [Wolfram Mathematica](http://www.wolfram.com/mathematica/) >=9 to read this performance data and generate plots from the preprint. Where plots are unannotated, Keynote was used to add back annotations.

All-of-GEUVADIS run
-----
1. Change the output bucket and argument of `--ec2-key-name` in `preprocess_all_of_GEUVADIS.sh`, and run the preprocess job flow contained in that script.
2. Change the output bucket and argument of `--ec2-key-name` in `submit_all_of_geuvadis_job.sh`, and run the alignment job flow contained in that script.

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
This generates `ssMat_geuvadis.rda` that has for all the 290416 ERs the percent of variance explained by each of the 16 variables: Population, RIN value, RNA extraction batch, RNA concentration, RNA quantity used, Library preparation date, Primer index, Method concentration measure, Library concentration, Library size, Library concentration used, Cluster kit, Sequencing kit, Cluster density, Lane, and Residual variation. It then creates boxplots saved in `r2_boxplots_overall.pdf` and a venn diagram of the ERs overlapping annotation features from ENSEMBL v75 saved in `ensemblVenn.pdf`. Overlap information is saved in `ensemblAnno.Rdata`. `derfinder_analysis.o4344170` and `derfinder_analysis.e4344170` are the corresponding log files.

GEUVADIS read count histogram (Figure 4 from preprint)
----
Download the 666 paired-end GEUVADIS sample FASTQ.gzs contained in `GEUVADIS_all_descriptive.manifest`. Run the command `gzip -cd | wc -l` on each file ending with `_1.fastq.gz` from the manifest, and divide the result by 2 to obtain the number of reads in the corresponding sample. We generated the histogram in [Wolfram Mathematica](http://www.wolfram.com/mathematica/) v10.
