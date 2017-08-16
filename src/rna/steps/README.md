* `preprocess.py` (optional)
    * Input: Manifest file
    * Output: `preproc`
* `align_reads.py`
    * Input: `preproc`
    * Output: `align_reads`
        * `align_reads/readletized`: Used by `align_readlets.py`
        * `align_reads/dummy`: Used by `junction_index.py`
        * `align_reads/unique`: Used by `cojunction_enum.py`
        * `align_reads/unmapped`: Used by `realign_reads.py`
        * `align_reads/sam`: Used by `compare_alignments.py`
        * `align_reads/postponed_sam`: Used by `bam.py`
        * `align_reads/exon_diff`: Used by `sum.py`
    * Calls: `align_reads_delegate.py`
    * Uses: Bowtie 2
        * `--sam-no-qname-trunc --local -t --no-hd --mm --12`
* `align_readlets.py`
    * Input: `align_reads/readletized`
    * Output: `align_readlets`
    * Calls: `align_readlets_delegate.py`
    * Uses: Bowtie
        * `-S -t --sam-nohead --mm --12`
* `junction_search.py`
    * Input: `align_readlets`
    * Output: `junction_search`
* `junction_filter.py`
    * Input: `junction_search`
    * Output: `junction_filter`
        * Output: `junction_filter/filter`: Used by `junction_config.py`
        * Output: `junction_filter/collect`: Used by `junction_collect.py`
* `junction_collect.py` (optional)
    * Input: `junction_filter/collect`
    * `--out`: `output/collected_junctions.tsv.gz`
* `junction_config.py`
    * Input: `junction_filter/filter`
    * Output: `junction_config`
* `junction_fasta.py`
    * Input: `junction_config`
    * Output: `junction_fasta`
* `junction_index.py`
    * Input: `junction_fasta`, `align_reads/dummy`
    * Output: `junction_index`
    * `--out`: `output/cross_sample_results`
* `cojunction_enum.py`
    * Input: `align_reads/unique`
    * Output: `cojunction_enum`
    * Calls: `cojunction_enum_delegate.py`
    * Uses: Bowtie 2
        * `--local -t --no-hd --mm --12 --score-min L,?,0 -D 24 -R 3 -N 1 -L 20 -i L,4,0`
* `cojunction_fasta.py`
    * Input: `cojunction_enum`
    * Output: `cojunction_fasta`
* `realign_reads.py`
    * Input: `align_reads/unmapped`, `cojunction_fasta`
    * Output: `realign_reads`
    * Calls: `realign_reads_delegate.py`
    * Uses: Bowtie 2
        * `-k ? --local -t --no-hd --mm -x`
* `compare_alignments.py`
    * Input: `align_reads/postponed_sam`, `realign_reads`
        * `compare_alignments/junction_bed`: Used by `junction_coverage.py`
        * `compare_alignments/sam_junction_ties`: Used by `junction_coverage.py`
        * `compare_alignments/sam_clip_ties`: Used by `break_ties.py`
        * `compare_alignments/sam`: Used by `bam.py`
        * `compare_alignments/exon_diff`: Used by `sum.py`
    * Output: `compare_alignments`
    * Uses: Bowtie??
* `junction_coverage.py`
    * Input: `compare_alignments/junction_bed`, `compare_alignments/sam_clip_ties`
    * Output: `junction_coverage`
* `break_ties.py`
    * Input: `junction_coverage`, `compare_alignments/sam_clip_ties`
    * Output: `junction_coverage`
        * `break_ties/sam`: Used by `bam.py`
        * `break_ties/exon_diff`: Used by `sum.py`
* `bam.py`
    * Input: `compare_alignments/sam`, `break_ties/sam`, `align_reads/sam`
    * Output: `bam`
        * `bam/counts`: Used by `collect_read_stats.py`
    * `--out`: `output/alignments`
* `collect_read_stats.py`
    * Input: `bam/counts`
    * Output: `read_counts`
    * `--out`: `output/cross_sample_results`
* `sum.py`
    * Input: `align_reads/exon_diff`, `compare_alignments/exon_diff`, `break_ties/exon_diff`
    * Output: `collapse`
* `coverage_pre.py`
    * Input: `collapse`
    * Files: `output/cross_sample_results/counts.tsv.gz`
    * Output: `precoverage`
* `coverage.py`
    * Input: `precoverage/coverage`
    * Output: `coverage`
    * `--out`: `output/coverage_bigwigs`
* `bed_pre.py`:
    * Input: `compare_alignments/indel_bed`, `break_ties/indel_bed`, `compare_alignments/junction_bed`, `break_ties/junction_bed`
    * Output: `prebed`
        * `prebed/collect`: Used by `tsv.py`
        * `prebed/bed`: Used by: `bed.py`
* `tsv.py`:
    * Input: `coverage`, `prebed/collect`
    * Output: `tsv`
* `bed.py`
    * Input: `prebed/bed`
    * Output: `bed`
