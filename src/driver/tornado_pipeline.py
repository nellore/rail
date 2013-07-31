"""
tornado_pipeline.py

"""

import pipeline

class PreprocessingStep(pipeline.Step):
    def __init__(self, inp, output, tconf, pconf):
        super(pipeline.Step, self).__init__(\
            "Preprocess", # name
            inp,      # input URL
            output,   # output URL
            "org.apache.hadoop.mapred.lib.NLineInputFormat",
            None,
            "python %s/preproc.py --compress=%s" % (pconf.emrLocalDir, tconf.preprocCompress))

class AlignStep(pipeline.Step):
    def __init__(self, inp, output, tconf, pconf):
        d = pconf.emrLocalDir
        super(pipeline.Step, self).__init__(\
            "Align", # name
            inp,     # input URL
            output,  # output URL
            None,    # input format
            None,    # no aggregation
            "python %s/align.py --refseq=%s/fasta/genome.fa --faidx=%s/fasta/genome.fa.fai --bowtieArgs='%s' --bowtieExe='%s' --bowtieIdx=%s/index/genome --readletLen %d --readletIval %d --partition-len %d" % (d, d, d, tconf.clusterRadius, tconf.bowtieArgs, tconf.bowtieExe, d, tconf.readletLen, tconf.readletIval, tconf.partitionLen))

class IntronStep(pipeline.Step):
    def __init__(self, inp, output, tconf, pconf):
        d = pconf.emrLocalDir
        super(pipeline.Step, self).__init__(\
            "Intron", # name
            inp,     # input URL
            output,  # output URL
            None,    # input format
            pipeline.Aggregation(None, 8, 1, 2),
            "cat",     # mapper
            "python %s/intron.py --refseq=%s/fasta/genome.fa --radius=%d --readletLen %d --readletIval %d" % (d, d, tconf.clusterRadius, tconf.readletLen, tconf.readletIval))

class MergeStep(pipeline.Step):
    def __init__(self, inp, output, tconf, pconf):
        d = pconf.emrLocalDir
        super(pipeline.Step, self).__init__(\
            "Merge", # name
            inp,     # input URL
            output,  # output URL
            None,    # input format
            pipeline.Aggregation(None, 8, 1, 2),
            "cat",     # mapper
            "python %s/merge.py" % d)

class WalkPrenormStep(pipeline.Step):
    def __init__(self, inp, output, tconf, pconf):
        d = pconf.emrLocalDir
        super(pipeline.Step, self).__init__(\
            "WalkPreNormalize", # name
            inp,     # input URL
            output,  # output URL
            None,    # input format
            pipeline.Aggregation(None, 8, 1, 2),
            "cat",     # mapper
            "python %s/walk_prenorm.py --partition-len %d" % (d, tconf.partitionLen))

class NormalizeStep(pipeline.Step):
    def __init__(self, inp, output, tconf, pconf):
        d = pconf.emrLocalDir
        super(pipeline.Step, self).__init__(\
            "Normalize", # name
            inp,     # input URL
            output,  # output URL
            None,    # input format
            pipeline.Aggregation(None, 4, 1, 1),
            "cat",     # mapper
            "python %s/normalize.py" % d)

class NormalizePostStep(pipeline.Step):
    def __init__(self, inp, output, tconf, pconf):
        super(pipeline.Step, self).__init__(\
            "NormalizePost", # name
            inp,     # input URL
            output,  # output URL
            None,    # input format
            pipeline.Aggregation(1, None, 0, 0),
            "cat",     # mapper
            "python %s/normalize_post.py" % emrLocalDir)

class WalkFitStep(pipeline.Step):
    def __init__(self, inp, output, tconf, pconf):
        super(pipeline.Step, self).__init__(\
            "WalkFit", # name
            inp,     # input URL
            output,  # output URL
            None,    # input format
            pipeline.Aggregation(None, 8, 1, 2),
            "cat",     # mapper
            "python %s/walk_fit.py" % emrLocalDir)

class EbayesStep(pipeline.Step):
    """ Just 1 reduce task """
    def __init__(self, inp, output, tconf, pconf):
        super(pipeline.Step, self).__init__(\
            "EmpiricalBayes", # name
            inp,     # input URL
            output,  # output URL
            None,    # input format
            pipeline.Aggregation(1, None, 0, 0),
            "cat",     # mapper
            "python %s/ebayes.py" % emrLocalDir)

class HmmParamsStep(pipeline.Step):
    def __init__(self, inp, output, tconf, pconf):
        super(pipeline.Step, self).__init__(\
            "HMMparams", # name
            inp,     # input URL
            output,  # output URL
            None,    # input format
            pipeline.Aggregation(1, None, 0, 0),
            "cat",     # mapper
            "python %s/hmm_params.py" % emrLocalDir)

class HmmStep(pipeline.Step):
    def __init__(self, inp, output, tconf, pconf):
        super(pipeline.Step, self).__init__(\
            "HMM", # name
            inp,     # input URL
            output,  # output URL
            None,    # input format
            pipeline.Aggregation(None, 8, 1, 2),
            "cat",     # mapper
            "python %s/hmm.py" % emrLocalDir)  # reducer

class AggrPathStep(pipeline.Step):
    def __init__(self, inp, output, tconf, pconf):
        super(pipeline.Step, self).__init__(\
            "HMMPaths", # name
            inp,     # input URL
            output,  # output URL
            None,    # input format
            pipeline.Aggregation(nperms * 2, None, 1, 2),
            "cat",     # mapper
            "python %s/aggr_path.py" % emrLocalDir)
