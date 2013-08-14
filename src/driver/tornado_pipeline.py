"""
tornado_pipeline.py

"""

import pipeline
import os
import site

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "util"))

import url

class PreprocessingStep(pipeline.Step):
    def __init__(self, inps, output, tconf, pconf):
        compressArg = ""
        if pconf.preprocCompress is not None and pconf.preprocCompress == "gzip":
            compressArg = "--gzip-output "
        super(PreprocessingStep, self).__init__(\
            inps,
            url.Url("hdfs:///dummy"),
            name="Preprocess", # name
            inputFormat="org.apache.hadoop.mapred.lib.NLineInputFormat",
            mapper="python %s/src/rnawesome/preprocess.py --nucs-per-file=10000000 %s--push=%s --ignore-first-token" % (pconf.emrLocalDir, compressArg, output.toUpperUrl()))

class AlignStep(pipeline.Step):
    def __init__(self, inps, output, tconf, pconf):
        d = pconf.emrLocalDir
        bexe = ""
        if tconf.bowtieExe is not None:
            bexe = "--bowtieExe='%s' " % tconf.bowtieExe
        super(AlignStep, self).__init__(\
            inps,
            output,
            name="Align",
            mapper="python %s/src/rnawesome/align.py --refseq=%s/fasta/genome.fa --faidx=%s/fasta/genome.fa.fai %s--bowtieIdx=%s/index/genome --readletLen %d --readletIval %d --partition-len %d -- %s" % (d, d, d, bexe, d, tconf.readletLen, tconf.readletIval, tconf.partitionLen, tconf.bowtieArgs()),
            outputFormat='edu.jhu.cs.MultipleOutputFormat',
            libjars=['/mnt/lib/multiplefiles.jar'])

class IntronStep(pipeline.Step):
    def __init__(self, inps, output, tconf, pconf):
        d = pconf.emrLocalDir
        super(IntronStep, self).__init__(\
            inps,
            output,  # output URL
            name="Intron", # name
            aggr=pipeline.Aggregation(None, 8, 1, 2), # 8 tasks per reducer
            reducer="python %s/src/rnawesome/intron.py --refseq=%s/fasta/genome.fa --radius=%d --readletLen=%d --readletIval=%d" % (d, d, tconf.clusterRadius, tconf.readletLen, tconf.readletIval))

class MergeStep(pipeline.Step):
    def __init__(self, inps, output, tconf, pconf):
        d = pconf.emrLocalDir
        super(MergeStep, self).__init__(\
            inps,
            output,  # output URL
            name="Merge", # name
            aggr=pipeline.Aggregation(None, 8, 1, 2),
            reducer="python %s/src/rnawesome/merge.py --partition-stats" % d,
            outputFormat='edu.jhu.cs.MultipleOutputFormat',
            libjars=['/mnt/lib/multiplefiles.jar'])

class WalkPrenormStep(pipeline.Step):
    def __init__(self, inps, output, tconf, pconf):
        d = pconf.emrLocalDir
        super(WalkPrenormStep, self).__init__(\
            inps,
            output,  # output URL
            name="WalkPreNormalize", # name
            aggr=pipeline.Aggregation(None, 8, 1, 2),
            reducer="python %s/src/rnawesome/walk_prenorm.py --partition-stats --partition-len=%d" % (d, tconf.partitionLen),
            outputFormat='edu.jhu.cs.MultipleOutputFormat',
            libjars=['/mnt/lib/multiplefiles.jar'])

class NormalizeStep(pipeline.Step):
    def __init__(self, inps, output, tconf, pconf):
        d = pconf.emrLocalDir
        super(NormalizeStep, self).__init__(\
            inps,
            output,  # output URL
            name="Normalize", # name
            aggr=pipeline.Aggregation(None, 4, 1, 1),
            reducer="python %s/src/rnawesome/normalize.py --percentile %f --out_dir=%s/coverage --bigbed_exe=%s/bin/bedToBigBed --faidx=%s/fasta/genome.fa.fai" % (d, tconf.normPercentile, pconf.out.toUpperUrl(), d, d))

class NormalizePostStep(pipeline.Step):
    def __init__(self, inps, output, tconf, pconf):
        d = pconf.emrLocalDir
        super(NormalizePostStep, self).__init__(\
            inps,
            output,  # output URL
            name="NormalizePost", # name
            aggr=pipeline.Aggregation(1, None, 0, 0),
            reducer="python %s/src/rnawesome/normalize_post.py --out=%s/normalize --manifest=%s/MANIFEST" % (d, pconf.out.toUpperUrl(), d))

class WalkFitStep(pipeline.Step):
    def __init__(self, inps, output, tconf, pconf):
        d = pconf.emrLocalDir
        super(WalkFitStep, self).__init__(\
            inps,
            output,  # output URL
            name="WalkFit", # name
            aggr=pipeline.Aggregation(None, 8, 1, 2),
            reducer="python %s/src/rnawesome/walk_fit.py" % d)

class EbayesStep(pipeline.Step):
    """ Just 1 reduce task """
    def __init__(self, inps, output, tconf, pconf):
        d = pconf.emrLocalDir
        super(EbayesStep, self).__init__(\
            inps,
            output,
            name="EmpiricalBayes", # name
            aggr=pipeline.Aggregation(1, None, 0, 0),
            reducer="python %s/src/rnawesome/ebayes.py" % d)

class HmmParamsStep(pipeline.Step):
    def __init__(self, inps, output, tconf, pconf):
        d = pconf.emrLocalDir
        super(HmmParamsStep, self).__init__(\
            inps,
            output,
            name="HMMparams", # name
            aggr=pipeline.Aggregation(1, None, 0, 0),
            reducer="python %s/src/rnawesome/hmm_params.py" % d)

class HmmStep(pipeline.Step):
    def __init__(self, inps, output, tconf, pconf):
        d = pconf.emrLocalDir
        super(HmmStep, self).__init__(\
            inps,
            output,
            name="HMM", # name
            aggr=pipeline.Aggregation(None, 8, 1, 2),
            reducer="python %s/src/rnawesome/hmm.py" % d)  # reducer

class AggrPathStep(pipeline.Step):
    def __init__(self, inps, output, tconf, pconf):
        d = pconf.emrLocalDir
        super(AggrPathStep, self).__init__(\
            inps,
            output,
            name="HMMPaths", # name
            aggr=pipeline.Aggregation(tconf.numPermutations * 2, None, 1, 2),
            reducer="python %s/src/rnawesome/aggr_path.py" % d)
