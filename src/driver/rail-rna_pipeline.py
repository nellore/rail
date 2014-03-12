"""
rail-rna_pipeline.py

"""

import pipeline
import os
import site
import re

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "util"))

import url
import version

#
# Note: All the %%(something)%% placeholders below are substituted in a
# mode-dependent manner.  That is, a different string could be substituted
# depending on whether we're in local, hadoop or emr mode.  This module tries
# to stay as mode-agnostic as possible.
#
# %%BASE%%: Base Rail directory
# %%MANIFEST%%: Manifest file
# %%BOWTIE%%: Path to Bowtie executable on workers/local machine
# %%BEDTOBIGBED%%: Path to bedToBigBed executable on workers/local machine
#

class PreprocessingStep(pipeline.Step):
    def __init__(self, inps, output, tconf, gconf):
        compressArg = ""
        if gconf.preprocCompress is not None and gconf.preprocCompress == "gzip":
            compressArg = "--gzip-output "
        mapperStr = """
            python %%BASE%%/src/rail-rna/preprocess.py 
                --nucs-per-file=8000000
                %s
                --push=%s
                --ignore-first-token""" % (compressArg, output.toUpperUrl())
        mapperStr = re.sub('\s+', ' ', mapperStr.strip())
        super(PreprocessingStep, self).__init__(\
            inps,
            None, # dummy: url.Url("hdfs:///dummy"),
            name="Preprocess", # name
            lineByLine=True,
            mapper=mapperStr)

class AlignStep(pipeline.Step):
    def __init__(self, inps, output, tconf, gconf):
        # For profiling, use the following mapperStr
        '''mapperStr = """
            python -m cProfile -o ~/align_output_profile %%BASE%%/src/rail-rna/align.py
                --refseq=%%REF_FASTA%% 
                --faidx=%%REF_FASTA_INDEX%% 
                --bowtie-idx=%%REF_BOWTIE_INDEX%% 
                --bowtie-exe=%%BOWTIE%% 
                --max-readlet-size %d 
                --readlet-interval %d 
                --partition-len %d 
                --exon-differentials 
                --verbose 
                -- %s""" % (tconf.readletLen, tconf.readletIval, tconf.partitionLen, tconf.bowtieArgs())'''
        mapperStr = """
            python %%BASE%%/src/rail-rna/align.py
                --bowtie-idx=%%REF_BOWTIE_INDEX%% 
                --bowtie-exe=%%BOWTIE%%
                --max-readlet-size %d 
                --readlet-interval %d 
                --partition-length %d
                --max-intron-size %d
                --min-intron-size %d
                --capping-fraction %f
                --exon-differentials 
                --verbose %s %s
                --min-cap-query-size %d
                --cap-search-window-size %d
                -- %s""" % (tconf.readletLen, tconf.readletIval, tconf.partitionLen,
                    tconf.max_intron_size,
                    tconf.min_intron_size, tconf.capping_fraction,
                    '--stranded' if tconf.stranded else '', 
                    '--do-not-search-for-caps' if tconf.do_not_search_for_caps else '',
                    tconf.min_cap_query_size,
                    tconf.cap_search_window_size,
                    tconf.bowtieAlignArgs())
        mapperStr = re.sub('\s+', ' ', mapperStr.strip())
        super(AlignStep, self).__init__(\
            inps,
            output,
            name="Align",
            mapper=mapperStr,
            multipleOutput=True)

class AlignPostStep(pipeline.Step):
    def __init__(self, inps, output, tconf, gconf):
        reducerStr = """
            python %%BASE%%/src/rail-rna/align_post.py
                --out=%s/spliced_alignments
                --refseq=%%REF_FASTA%%
                --samtools-exe=%%SAMTOOLS%%
                --bam-basename=%s
                %s
                %s
            """ % (gconf.out, tconf.bam_basename, 
                    '--output-by-chromosome' if tconf.output_bam_by_chromosome else '',
                    '--output-sam' if tconf.output_sam else '')
        reducerStr = re.sub('\s+', ' ', reducerStr.strip())
        super(AlignPostStep, self).__init__(\
            inps,
            output,
            name="AlignPost",
            aggr=(pipeline.Aggregation(None, 1, 1, 2) if (tconf.output_bam_by_chromosome and gconf.out is not None) \
                    else pipeline.Aggregation(1, None, 1, 2)),
            reducer=reducerStr)

class IntronStep(pipeline.Step):
    def __init__(self, inps, output, tconf, gconf):
        reducerStr = """
            python %%BASE%%/src/rail-rna/intron.py
                --bowtie-idx=%%REF_BOWTIE_INDEX%% 
                --cluster-radius=%d
                --intron-partition-overlap=%d
                --verbose
                --partition-length %d
                --motif-radius=%d
                %s
        """ % (tconf.clusterRadius, tconf.intronPartitionOlap,
                tconf.partitionLen, tconf.motifRadius,
                '--stranded' if tconf.stranded else '')
        reducerStr = re.sub('\s+', ' ', reducerStr.strip())
        super(IntronStep, self).__init__(\
            inps,
            output,  # output URL
            name="Intron",  # name
            aggr=pipeline.Aggregation(None, 8, 1, 2),  # 8 tasks per reducer
            reducer=reducerStr,
            multipleOutput=True)

class IntronPostStep(pipeline.Step):
    def __init__(self, inps, output, tconf, gconf):
        reducerStr = """
            python %%BASE%%/src/rail-rna/intron_post.py
                --bowtie-idx=%%REF_BOWTIE_INDEX%% 
                --out=%s/index
            """ % gconf.intermediate
        reducerStr = re.sub('\s+', ' ', reducerStr.strip())
        super(IntronPostStep, self).__init__(\
            inps,
            output,
            name="IntronPost",
            aggr=pipeline.Aggregation(1, None, 1, 5),
            reducer=reducerStr)

class MergeStep(pipeline.Step):
    def __init__(self, inps, output, tconf, gconf):
        reducerStr = "python %%BASE%%/src/rail-rna/merge.py --partition-stats"
        super(MergeStep, self).__init__(\
            inps,
            output,  # output URL
            name="Merge", # name
            aggr=pipeline.Aggregation(None, 8, 1, 2),
            reducer=reducerStr,
            multipleOutput=True)

class RealignStep(pipeline.Step):
    def __init__(self, inps, output, tconf, gconf):
        mapperStr = """
            python %%BASE%%/src/rail-rna/realign.py
                --bowtie-idx=%s/index/intron
                --bowtie-exe=%%BOWTIE%%
                --partition-length %d
                --exon-differentials 
                --verbose %s
                -- %s""" % (gconf.intermediate, tconf.partitionLen,
                    '--stranded' if tconf.stranded else '',
                    tconf.bowtieRealignArgs())
        mapperStr = re.sub('\s+', ' ', mapperStr.strip())
        super(RealignStep, self).__init__(\
            inps,
            output,
            name="Realign",
            mapper=mapperStr,
            multipleOutput=True)

class WalkPrenormStep(pipeline.Step):
    def __init__(self, inps, output, tconf, gconf):
        reducerStr = """
            python %%BASE%%/src/rail-rna/walk_prenorm.py 
                --partition-stats 
                --partition-len=%d""" % (tconf.partitionLen)
        reducerStr = re.sub('\s+', ' ', reducerStr.strip())
        super(WalkPrenormStep, self).__init__(\
            inps,
            output,  # output URL
            name="WalkPreNormalize", # name
            aggr=pipeline.Aggregation(None, 8, 1, 2),
            reducer=reducerStr,
            multipleOutput=True)

class NormalizePreStep(pipeline.Step):
    def __init__(self, inps, output, tconf, gconf):
        reducerStr = """
            python %%BASE%%/src/rail-rna/normalize_pre.py 
                --partition-stats 
                --partition-len=%d""" % (tconf.partitionLen)
        reducerStr = re.sub('\s+', ' ', reducerStr.strip())
        super(NormalizePreStep, self).__init__(\
            inps,
            output,  # output URL
            name="NormalizePre", # name
            aggr=pipeline.Aggregation(None, 8, 2, 3),
            reducer=reducerStr,
            multipleOutput=True)

class CollapseStep(pipeline.Step):
    def __init__(self, inps, output, tconf, gconf):
        reducerStr = """
            python %%BASE%%/src/rail-rna/collapse.py 
                %s
            """ % ('--stranded' if tconf.stranded else '')
        reducerStr = re.sub('\s+', ' ', reducerStr.strip())
        super(CollapseStep, self).__init__(\
            inps,
            output,  # output URL
            name="Collapse", # name
            aggr=pipeline.Aggregation(None, 8, 1, 1),
            reducer=reducerStr,
            multipleOutput=True)

class CoveragePreStep(pipeline.Step):
    def __init__(self, inps, output, tconf, gconf):
        reducerStr = """
            python %%BASE%%/src/rail-rna/coverage_pre.py 
                --partition-stats 
                --partition-len=%d""" % (tconf.partitionLen)
        reducerStr = re.sub('\s+', ' ', reducerStr.strip())
        super(CoveragePreStep, self).__init__(\
            inps,
            output,  # output URL
            name="CoveragePre", # name
            aggr=pipeline.Aggregation(None, 8, 2, 3),
            reducer=reducerStr,
            multipleOutput=True)

class CoverageStep(pipeline.Step):
    def __init__(self, inps, output, tconf, gconf):
        reducerStr = """
            python %%BASE%%/src/rail-rna/coverage.py 
                --percentile %f
                --out=%s/coverage
                --bigbed-exe=%%BEDTOBIGBED%%
                --refseq=%%REF_FASTA%% 
                --faidx=%%REF_FASTA_INDEX%%
                --verbose""" % (tconf.normPercentile, gconf.out)
        reducerStr = re.sub('\s+', ' ', reducerStr.strip())
        super(CoverageStep, self).__init__(\
            inps,
            output,  # output URL
            name="Coverage", # name
            aggr=pipeline.Aggregation(None, 1, 1, 3),
            reducer=reducerStr)

class CoveragePostStep(pipeline.Step):
    def __init__(self, inps, output, tconf, gconf):
        reducerStr = """
            python %%BASE%%/src/rail-rna/coverage_post.py 
                --out=%s/normalize 
                --manifest=%%MANIFEST%%""" % gconf.out
        reducerStr = re.sub('\s+', ' ', reducerStr.strip())
        super(CoveragePostStep, self).__init__(\
            inps,
            output,  # output URL
            name="CoveragePost", # name
            aggr=pipeline.Aggregation(1, None, 0, 0),
            reducer=reducerStr)

class NormalizeStep(pipeline.Step):
    def __init__(self, inps, output, tconf, gconf):
        reducerStr = """
            python %%BASE%%/src/rail-rna/normalize2.py 
                --percentile %f
                --out_dir=%s/coverage
                --bigbed_exe=%%BEDTOBIGBED%%
                --faidx=%%REF_FASTA_INDEX%%
                --verbose""" % (tconf.normPercentile, gconf.out)
        reducerStr = re.sub('\s+', ' ', reducerStr.strip())
        super(NormalizeStep, self).__init__(\
            inps,
            output,  # output URL
            name="Normalize", # name
            aggr=pipeline.Aggregation(None, 1, 1, 3),
            reducer=reducerStr)

class NormalizePostStep(pipeline.Step):
    def __init__(self, inps, output, tconf, gconf):
        reducerStr = """
            python %%BASE%%/src/rail-rna/normalize_post.py 
                --out=%s/normalize 
                --manifest=%%MANIFEST%%""" % gconf.out
        reducerStr = re.sub('\s+', ' ', reducerStr.strip())
        super(NormalizePostStep, self).__init__(\
            inps,
            output,  # output URL
            name="NormalizePost", # name
            aggr=pipeline.Aggregation(1, None, 0, 0),
            reducer=reducerStr)

class WalkFitStep(pipeline.Step):
    def __init__(self, inps, output, tconf, gconf):
        reducerStr = """
            python %%BASE%%/src/rail-rna/walk_fit.py"""
        reducerStr = re.sub('\s+', ' ', reducerStr.strip())
        super(WalkFitStep, self).__init__(\
            inps,
            output,  # output URL
            name="WalkFit", # name
            aggr=pipeline.Aggregation(None, 8, 1, 2),
            reducer=reducerStr)

class EbayesStep(pipeline.Step):
    """ Just 1 reduce task """
    def __init__(self, inps, output, tconf, gconf):
        reducerStr = """
            python %%BASE%%/src/rail-rna/walk_fit.py"""
        reducerStr = re.sub('\s+', ' ', reducerStr.strip())
        super(EbayesStep, self).__init__(\
            inps,
            output,
            name="EmpiricalBayes", # name
            aggr=pipeline.Aggregation(1, None, 0, 0),
            reducer=reducerStr)

class HmmParamsStep(pipeline.Step):
    def __init__(self, inps, output, tconf, gconf):
        reducerStr = """
            python %%BASE%%/src/rail-rna/hmm_params.py"""
        reducerStr = re.sub('\s+', ' ', reducerStr.strip())
        super(HmmParamsStep, self).__init__(\
            inps,
            output,
            name="HMMparams", # name
            aggr=pipeline.Aggregation(1, None, 0, 0),
            reducer=reducerStr)

class HmmStep(pipeline.Step):
    def __init__(self, inps, output, tconf, gconf):
        reducerStr = """
            python %%BASE%%/src/rail-rna/hmm.py"""
        reducerStr = re.sub('\s+', ' ', reducerStr.strip())
        super(HmmStep, self).__init__(\
            inps,
            output,
            name="HMM", # name
            aggr=pipeline.Aggregation(None, 8, 1, 2),
            reducer=reducerStr)

class AggrPathStep(pipeline.Step):
    def __init__(self, inps, output, tconf, gconf):
        reducerStr = """
            python %%BASE%%/src/rail-rna/aggr_path.py"""
        reducerStr = re.sub('\s+', ' ', reducerStr.strip())
        super(AggrPathStep, self).__init__(\
            inps,
            output,
            name="HMMPaths", # name
            aggr=pipeline.Aggregation(tconf.numPermutations * 2, None, 1, 2),
            reducer=reducerStr)
