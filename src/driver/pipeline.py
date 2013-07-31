"""
pipeline.py

Classes and routines for describing steps in the pipeline and how they can be
implemented as MapReduce steps.
"""

class Aggregation(object):
    """ Encapsulates information about the aggregation before a reducer. """
    
    def __init__(self, ntasks, ntasksPerReducer, nbin, nsort):
        self.ntasks = ntasks
        self.ntasksPerReducer = ntasksPerReducer
        assert nsort >= nbin
        self.nbin = nbin
        self.nsort = nsort
    
    def toHadoopArgs(self, config):
        begArgs, endArgs = [], []
        v = config.hadoopVersion
        confArg = "-D"
        if v[0] == 0 and v[1] < 19:
            confArg = "-jobconf"
        if self.ntasks is not None:
            begArgs.append('"%s", "mapred.reduce.tasks=%d",' % (confArg, self.ntasks))
        else:
            assert self.ntasksPerReducer is not None
            begArgs.append('"%s", "mapred.reduce.tasks=%d",' % (confArg, self.ntasksPerReducer * config.nReducers))
            if v[0] == 0 and v[1] < 19:
                begArgs.append('"%s", "num.key.fields.for.partition=%d",' % (confArg, self.nbin))
            else:
                begArgs.append('"%s", "mapred.text.key.partitioner.options=-k1,%d",' % (confArg, self.nbin))
            begArgs.append('"%s", "stream.num.map.output.key.fields=%d",' % (confArg, self.nsort))
            endArgs.append('"-partitioner", "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",')
        return begArgs, endArgs

class PipelineConfig(object):
    """ Information about the particular system we're going to run on, which
        changes some of the exact arguments we use """
    
    def __init__(self, hadoopVersion, waitOnFail, emrStreamJar, nReducers, emrLocalDir, preprocCompress):
        # Hadoop version used on EMR cluster
        self.hadoopVersion = hadoopVersion
        # Whether to keep the EMR cluster running if job fails
        self.waitOnFail = waitOnFail
        # Path to streaming jar
        self.emrStreamJar = emrStreamJar
        # Total number of reducers that can run simultaneously
        self.nReducers = nReducers
        # Local directory where reference jar and scripts have been installed
        self.emrLocalDir = emrLocalDir
        # Whether & how to compress output from preprocessor
        self.preprocCompress = preprocCompress

class Step(object):
    """ Encapsulates a single step of the pipeline, i.e. a single
        MapReduce/Hadoop job """
    
    def __init__(self, name, inp, output, inputFormat, aggr, mapper, reducer, out):
        self.name = name
        self.input = inp
        self.output = output
        self.inputFormat = inputFormat
        self.aggr = aggr
        self.mapper = mapper
        self.reducer = reducer
        self.out = out
    
    def toHadoopCmd(self):
        raise RuntimeError("toHadoopCmd not yet implemented")
    
    def toEmrCmd(self, config):
        lines = []
        lines.append('{')
        lines.append('  "Name" : "%s",' % self.name)
        lines.append('  "ActionOnFailure" : "%s",' % ("CANCEL_AND_WAIT" if config.waitOnFail else "TERMINATE_JOB_FLOW"))
        lines.append('  "HadoopJarStep": {')
        lines.append('    "Jar": "%s",' % config.emrStreamJar)
        lines.append('    "Args": [')
        
        begArgs, endArgs = self.reduceStyle.toHadoopArgs(config)
        for a in begArgs:
            lines.append('      ' + a)
        
        lines.append('      "-input", "%s",' % self.input)
        lines.append('      "-output", "%s",' % self.output)
        lines.append('      "-mapper", "%s",' % self.mapper)
        
        if self.reducer is not None:
            lines.append('      "-reducer", "%s",' % self.reducer)
        
        if self.inputFormat is not None:
            lines.append('      "-inputformat", "%s",' % self.inputFormat)
        
        for a in endArgs:
            lines.append('      ' + a)
        
        # Add cache files for this job
        for fn in self.cacheFiles:
            lines.append('      "-file", "%s",' % fn)
        
        # Add cache files for any tools used by this job
        
        lines.append('    ]')
        lines.append('  }')
        lines.append('}')
        
        return '\n'.join(lines)

