#!/usr/bin/env bash
## A script that:
## (1) Compresses the FASTQs and stages them on S3
## (2) Creates a manifest file for the FASTQs
## (3) Runs Rail-RNA on the 112 datasets on Elastic MapReduce WITH AND WITHOUT intron filter of 0.5,5
## $1: directory with fastqs
## $2: where on S3 to stage fastq.gzs
## $3: where on S3 to output Rail results
## $4: path to manifest file to write (taken in our experiments to be eval/GEUVADIS_112_sim.manifest)
## Requires AWS CLI
## Full command we ran was sh run_112_sim.sh /scratch2/langmead-fs1/geuvadis_sims_for_paper_v2 s3://rail-eu-west-1/geuv112sim_v2 s3://rail-eu-west-1/geuv112sim_v2.out /scratch0/langmead-fs1/rail/eval/GEUVADIS_112_sim.manifest
## Make sure Rail-RNA v0.1.9 is used to reproduce our results
FASTQDIR=$1
S3STAGED=$2
S3DEST=$3
MANIFEST=$4
cd $FASTQDIR
for fastq in $(find . -name \*.fastq ! -name \*right\*.fastq ! -name \*left\*.fastq | cut -c3-)
do
	leftfile=$(echo ${fastq} | rev | cut -c 7- | rev)_left.fastq
	rightfile=$(echo ${fastq} | rev | cut -c 7- | rev)_right.fastq
	awk '(NR-1) % 8 < 4' $fastq >$leftfile
	awk '(NR-1) % 8 >= 4' $fastq >$rightfile
done
find . -name \*.fastq ! -name \*right\*.fastq ! -name \*left\*.fastq | cut -c3- | python -c "import sys
for line in sys.stdin:
    sample_name = list(line.strip().rpartition('_')[0].partition('-'))
    sample_name[0] = sample_name[0] + '_sim'
    sample_root = line.strip().rpartition('.')[0]
    print '\t'.join(['${S3STAGED}/' + sample_root + '_left.fastq', '0', '${S3STAGED}/' + sample_root + '_right.fastq', '0', ''.join(sample_name)])" >$MANIFEST
for i in $(cat <(find . -name \*left.fastq) <(find . -name \*right.fastq) | awk '{print substr($0, 3)}'); do aws s3 cp $i $S3STAGED/$i; done
echo "Now execute these commands to submit the jobs."
echo "rail-rna go elastic -c 40 -a hg19 -m $MANIFEST -o $S3DEST --core-instance-bid-price 0.11 --master-instance-bid-price 0.11 --region eu-west-1"
# WITHOUT intron filter
echo "rail-rna go elastic -c 40 -a hg19 -m $MANIFEST -o $S3DEST.nofilter --core-instance-bid-price 0.11 --master-instance-bid-price 0.11 --region eu-west-1 --intron-criteria 0,0"
