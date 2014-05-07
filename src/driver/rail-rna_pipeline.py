"""
rail-rna_pipeline.py
"""

import pipeline
import os
import site
import re

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "util"))

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
    def __init__(self, inps, output, _, gconf):
        compress_arg = ""
        if gconf.preprocCompress is not None and gconf.preprocCompress == "gzip":
            compress_arg = "--gzip-output "
        push_output = output.plus('push')
        mapper_str = """
            %%PYPY%% %%BASE%%/src/rail-rna/preprocess.py 
                --nucs-per-file=8000000
                %s
                --push=%s
                --ignore-first-token""" % (compress_arg, push_output.toUpperUrl())
        mapper_str = re.sub('\s+', ' ', mapper_str.strip())
        super(PreprocessingStep, self).__init__(
            inps,
            output,
            name="Preprocess",  # name
            lineByLine=True,
            mapper=mapper_str)


class AlignReadsStep(pipeline.Step):
    def __init__(self, inps, output, tconf, _):
        mapper_str = """
            %%PYPY%% %%BASE%%/src/rail-rna/align_reads.py
                --bowtie-idx=%%REF_BOWTIE_INDEX%%
                --bowtie2-idx=%%REF_BOWTIE2_INDEX%%
                --bowtie2-exe=%%BOWTIE2%%
                --exon-differentials
                --partition-length %d
                --manifest=%%MANIFEST%%
                --verbose
                -- %s --local """ % (tconf.partitionLen,
                             tconf.bowtieArgs())
        combiner_str = """
            %%PYPY%% %%BASE%%/src/rail-rna/sum.py
                --type 3
                --value-count 2
                --target %s
        """ % ('readletize',)
        mapper_str = re.sub('\s+', ' ', mapper_str.strip())
        combiner_str = re.sub('\s+', ' ', combiner_str.strip())
        super(AlignReadsStep, self).__init__(
            inps,
            output,
            name="AlignReads",
            mapper=mapper_str,
            combiner=combiner_str,
            multipleOutput=True)

class CombineSequencesStep(pipeline.Step):
    def __init__(self, inps, output, tconf, _):
        reducer_str = """
            %%PYPY%% %%BASE%%/src/rail-rna/sum.py
                --type 3
                --value-count 2
                %s
        """ % ('',)
        reducer_str = re.sub('\s+', ' ', reducer_str.strip())
        super(CombineSequencesStep, self).__init__(
            inps,
            output,  # output URL
            name="CombineSequences",  # name
            aggr=pipeline.Aggregation(None, 4, 1, 1),  # 4 tasks per reducer
            reducer=reducer_str)

class ReadletizeStep(pipeline.Step):
    def __init__(self, inps, output, tconf, _):
        reducer_str = """
            %%PYPY%% %%BASE%%/src/rail-rna/readletize.py
                --max-readlet-size %d 
                --readlet-interval %d 
                --capping-multiplier %f
        """ % (tconf.readletLen, tconf.readletIval, tconf.capping_multiplier)
        combiner_str = """
            %%PYPY%% %%BASE%%/src/rail-rna/sum.py
                --type 3 %s
        """ % ('',)
        reducer_str = re.sub('\s+', ' ', reducer_str.strip())
        combiner_str = re.sub('\s+', ' ', combiner_str.strip())
        super(ReadletizeStep, self).__init__(
            inps,
            output,  # output URL
            name="Readletize",  # name
            aggr=pipeline.Aggregation(None, 4, 1, 1),  # 4 tasks per reducer
            reducer=reducer_str,
            combiner=combiner_str)

class CombineSubsequencesStep(pipeline.Step):
    def __init__(self, inps, output, tconf, _):
        reducer_str = """
            %%PYPY%% %%BASE%%/src/rail-rna/sum.py
                --type 3
                %s
        """ % ('',)
        reducer_str = re.sub('\s+', ' ', reducer_str.strip())
        super(CombineSubsequencesStep, self).__init__(
            inps,
            output,  # output URL
            name="CombineSubsequences",  # name
            aggr=pipeline.Aggregation(None, 4, 1, 1),  # 4 tasks per reducer
            reducer=reducer_str)

class AlignReadletsStep(pipeline.Step):
    def __init__(self, inps, output, tconf, _):
        reducer_str = """
            %%PYPY%% %%BASE%%/src/rail-rna/align_readlets.py
                --bowtie-idx=%%REF_BOWTIE_INDEX%% 
                --bowtie-exe=%%BOWTIE%%
                --verbose
                -- %s
        """ % ('-t --sam-nohead --startverbose -v 0 -a -m 80',)
        reducer_str = re.sub('\s+', ' ', reducer_str.strip())
        super(AlignReadletsStep, self).__init__(
            inps,
            output,  # output URL
            name="AlignReadlets",  # name
            aggr=pipeline.Aggregation(None, 4, 1, 1),  # 4 tasks per reducer
            reducer=reducer_str)

class IntronSearchStep(pipeline.Step):
    def __init__(self, inps, output, tconf, _):
        reducer_str = """
            %%PYPY%% %%BASE%%/src/rail-rna/intron_search.py
                --bowtie-idx=%%REF_BOWTIE_INDEX%% 
                --partition-length %d
                --max-intron-size %d
                --min-intron-size %d
                --min-exon-size %d
                --search-window-size %d
                --motif-radius %d
                --verbose
                """ % (tconf.partitionLen,
                            tconf.max_intron_size,
                            tconf.min_intron_size,
                            tconf.min_exon_size,
                            tconf.search_window_size,
                            tconf.motif_radius
                        )
        reducer_str = re.sub('\s+', ' ', reducer_str.strip())
        super(IntronSearchStep, self).__init__(
            inps,
            output,  # output URL
            name="IntronSearch",  # name
            aggr=pipeline.Aggregation(None, 4, 1, 1),  # 4 tasks per reducer
            reducer=reducer_str)

class IntronConfigStep(pipeline.Step):
    def __init__(self, inps, output, tconf, gconf):
        reducer_str = """
            %%PYPY%% %%BASE%%/src/rail-rna/intron_config.py
                --readlet-size %d
                --verbose
            """ % (tconf.readletLen,)
        reducer_str = re.sub('\s+', ' ', reducer_str.strip())
        super(IntronConfigStep, self).__init__(
            inps,
            output,
            name="IntronConfig",
            aggr=pipeline.Aggregation(None, 1, 2, 4),
            reducer=reducer_str)

class IntronFastaStep(pipeline.Step):
    def __init__(self, inps, output, tconf, gconf):
        reducer_str = """
            %%PYPY%% %%BASE%%/src/rail-rna/intron_fasta.py
                --verbose
                --bowtie-idx=%%REF_BOWTIE_INDEX%%
                %s
            """ % ('',)
        reducer_str = re.sub('\s+', ' ', reducer_str.strip())
        super(IntronFastaStep, self).__init__(
            inps,
            output,
            name="IntronFasta",
            aggr=pipeline.Aggregation(None, 8, 4, 4),
            reducer=reducer_str)

class IntronIndexStep(pipeline.Step):
    def __init__(self, inps, output, tconf, gconf):
        reducer_str = """
            %%PYPY%% %%BASE%%/src/rail-rna/intron_index.py
                --bowtie-build-exe=%%BOWTIE-BUILD%%
                --bowtie-idx=%%REF_BOWTIE_INDEX%%
                --out=%s/index
                %s
            """ % (gconf.out, '--keep-alive' if tconf.keep_alive else '')
        reducer_str = re.sub('\s+', ' ', reducer_str.strip())
        super(IntronIndexStep, self).__init__(
            inps,
            output,
            name="IntronIndex",
            aggr=pipeline.Aggregation(1, None, 1, 1),
            reducer=reducer_str)

class RealignReadletsStep(pipeline.Step):
    def __init__(self, inps, output, tconf, gconf, cache=None):
        reducer_str = """
            %%PYPY%% %%BASE%%/src/rail-rna/align_readlets.py
                --bowtie-idx=%%REF_INTRON_INDEX%%
                --bowtie-exe=%%BOWTIE%%
                --verbose
                -- %s
        """ % ('-t --sam-nohead --startverbose -v 0 -a -m 80',)
        reducer_str = re.sub('\s+', ' ', reducer_str.strip())
        super(RealignReadletsStep, self).__init__(
            inps,
            output,  # output URL
            name="RealignReadlets",  # name
            aggr=pipeline.Aggregation(None, 4, 1, 1),  # 4 tasks per reducer
            reducer=reducer_str,
            cache=cache)

class CointronSearchStep(pipeline.Step):
    def __init__(self, inps, output, tconf, _):
        reducer_str = """
            %%PYPY%% %%BASE%%/src/rail-rna/cointron_search.py
                --verbose %s
                """ % ('',)
        reducer_str = re.sub('\s+', ' ', reducer_str.strip())
        super(CointronSearchStep, self).__init__(
            inps,
            output,  # output URL
            name="CointronSearch",  # name
            aggr=pipeline.Aggregation(None, 4, 1, 1),  # 4 tasks per reducer
            reducer=reducer_str)

class CointronFastaStep(pipeline.Step):
    def __init__(self, inps, output, tconf, gconf):
        reducer_str = """
            %%PYPY%% %%BASE%%/src/rail-rna/cointron_fasta.py
                --verbose
                --bowtie-idx=%%REF_BOWTIE_INDEX%%
                %s
            """ % ('',)
        reducer_str = re.sub('\s+', ' ', reducer_str.strip())
        super(CointronFastaStep, self).__init__(
            inps,
            output,
            name="CointronFasta",
            aggr=pipeline.Aggregation(None, 8, 2, 2),
            reducer=reducer_str)

class RealignReadsStep(pipeline.Step):
    def __init__(self, inps, output, tconf, gconf, cache=None):
        reducer_str = """
            %%PYPY%% %%BASE%%/src/rail-rna/realign_reads.py
                --original-idx=%%REF_BOWTIE_INDEX%%
                --bowtie2-exe=%%BOWTIE2%%
                --partition-length %d
                --exon-differentials
                --manifest=%%MANIFEST%%
                --verbose %s %s
                -- %s --end-to-end """ % (tconf.partitionLen,
                            '--stranded' if tconf.stranded else '',
                            '--keep-alive' if tconf.keep_alive else '',
                            tconf.bowtieArgs())
        combiner_str = """
            %%PYPY%% %%BASE%%/src/rail-rna/sum.py 
                --target %s
            """ % ('exon_diff',)
        reducer_str = re.sub('\s+', ' ', reducer_str.strip())
        combiner_str = re.sub('\s+', ' ', combiner_str.strip())
        super(RealignReadsStep, self).__init__(
            inps,
            output,
            name="RealignReads",
            aggr=pipeline.Aggregation(None, 4, 1, 1),  # 4 tasks per reducer
            reducer=reducer_str,
            combiner=combiner_str,
            multipleOutput=True,
            )

class CollapseStep(pipeline.Step):
    def __init__(self, inps, output, tconf, _):
        reducer_str = """
            %%PYPY%% %%BASE%%/src/rail-rna/sum.py 
                %s
            """ % ('',)
        reducer_str = re.sub('\s+', ' ', reducer_str.strip())
        super(CollapseStep, self).__init__(
            inps,
            output,  # output URL
            name="Collapse",  # name
            aggr=pipeline.Aggregation(None, 8, 1, 1),
            reducer=reducer_str)


class CoveragePreStep(pipeline.Step):
    def __init__(self, inps, output, tconf, _):
        reducer_str = """
            %%PYPY%% %%BASE%%/src/rail-rna/coverage_pre.py
                --bowtie-idx=%%REF_BOWTIE_INDEX%%
                --partition-stats %s""" % ('',)
        reducer_str = re.sub('\s+', ' ', reducer_str.strip())
        super(CoveragePreStep, self).__init__(
            inps,
            output,  # output URL
            name="CoveragePre",  # name
            aggr=pipeline.Aggregation(None, 8, 2, 3),
            reducer=reducer_str,
            multipleOutput=True)


class CoverageStep(pipeline.Step):
    def __init__(self, inps, output, tconf, gconf):
        reducer_str = """
            %%PYPY%% %%BASE%%/src/rail-rna/coverage.py
                --bowtie-idx=%%REF_BOWTIE_INDEX%%
                --percentile %f
                --out=%s/coverage
                --bigbed-exe=%%BEDTOBIGBED%%
                --manifest=%%MANIFEST%%
                %s
                --verbose """ % (tconf.normPercentile, gconf.out,
                                '--keep-alive' if tconf.keep_alive else '',
                                )
        reducer_str = re.sub('\s+', ' ', reducer_str.strip())
        super(CoverageStep, self).__init__(
            inps,
            output,  # output URL
            name="Coverage",  # name
            aggr=pipeline.Aggregation(None, 1, 1, 3),
            reducer=reducer_str)


class CoveragePostStep(pipeline.Step):
    def __init__(self, inps, output, _, gconf):
        reducer_str = """
            %%PYPY%% %%BASE%%/src/rail-rna/coverage_post.py 
                --out=%s/normalize 
                --manifest=%%MANIFEST%%""" % gconf.out
        reducer_str = re.sub('\s+', ' ', reducer_str.strip())
        super(CoveragePostStep, self).__init__(
            inps,
            output,  # output URL
            name="CoveragePost",  # name
            aggr=pipeline.Aggregation(1, None, 1, 2),
            reducer=reducer_str)

class BamStep(pipeline.Step):
    def __init__(self, inps, output, tconf, gconf):
        reducer_str = """
            %%PYPY%% %%BASE%%/src/rail-rna/bam.py
                --out=%s/bam
                --bowtie-idx=%%REF_BOWTIE_INDEX%% 
                --samtools-exe=%%SAMTOOLS%%
                --bam-basename=%s
                --manifest=%%MANIFEST%%
                %s
                %s
                %s
            """ % (gconf.out, tconf.bam_basename,
                   '--output-by-chromosome' if tconf.output_bam_by_chromosome else '',
                   '--output-sam' if tconf.output_sam else '',
                   '--keep-alive' if tconf.keep_alive else '')
        reducer_str = re.sub('\s+', ' ', reducer_str.strip())
        super(BamStep, self).__init__(
            inps,
            output,
            name="Bam",
            aggr=(pipeline.Aggregation(None, 1, (2 if tconf.output_bam_by_chromosome else 1), 3)
                  if gconf.out is not None else pipeline.Aggregation(None, 1, 1, 3)),
            reducer=reducer_str)


class BedPreStep(pipeline.Step):
    def __init__(self, inps, output, _, _2):
        reducer_str = """
            %%PYPY%% %%BASE%%/src/rail-rna/bed_pre.py
                %s
        """ % ''
        reducer_str = re.sub('\s+', ' ', reducer_str.strip())
        super(BedPreStep, self).__init__(
            inps,
            output,  # output URL
            name="BedPre",  # name
            aggr=pipeline.Aggregation(None, 8, 6, 6),  # 8 tasks per reducer
            reducer=reducer_str)


class BedStep(pipeline.Step):
    def __init__(self, inps, output, tconf, gconf):
        reducer_str = """
            %%PYPY%% %%BASE%%/src/rail-rna/bed.py 
                --bowtie-idx=%%REF_BOWTIE_INDEX%%
                --out=%s/bed
                --manifest=%%MANIFEST%%
                --bed-basename=%s""" % (gconf.out, tconf.bed_basename)
        reducer_str = re.sub('\s+', ' ', reducer_str.strip())
        super(BedStep, self).__init__(
            inps,
            output,  # output URL
            name="Bed",  # namef
            aggr=pipeline.Aggregation(None, 1, 2, 4),
            reducer=reducer_str)