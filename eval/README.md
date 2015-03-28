Reproducing preprint results
-----
To reproduce results from the [preprint](http://finishit.com), perform the following steps. Note that input and output directories in scripts may need to be changed.
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
with 56 datasets, 8 from each lab, and a random sample C of B with 28 datasets, 4 from each lab. (Again, see script for data sources.) The corresponding manifest files are stored in `GEUVADIS_112.manifest`, `GEUVADIS_56.manifest`, and `GEUVADIS_28.manifest`. Optionally, use `make_sample_labels_consistent.py` to ensure that sample labels from all manifest files are the same for the same datasets. (Otherwise, replicate numbers could differ between two manifests for the same dataset).
3. Change the output bucket and home directory of Rail-RNA's source code (`rail/src`) in `prep_scaling_experiments.sh`. Run
```
sh prep_scaling_experiments.sh
```
to submit all jobs that preprocess GEUVADIS data for scaling experiments to Elastic MapReduce. This requires an account with [Amazon Web Services](http://aws.amazon.com/) and costs money. Note also that `--ec2-key-name` in all commands from `prep_scaling_experiments.sh` should either be removed to changed to a valid EC2 key name. The key name permits SSHing to EC2 instances while the preprocessing jobs are running.
4. 
Accuracy experiments
-----
Prerequisite: manifest file GEUVADIS_112.manifest from the scaling experiments above.