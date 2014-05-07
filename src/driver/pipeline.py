"""
pipeline.py

Classes and routines for describing steps in the pipeline and how they can be
implemented as MapReduce steps.

A few pieces of vocabulary:
- The overall DAG is a "flow"
- Each straight-line piece of the DAG is a "pipeline"
- Each Map-Aggregate-Reduce is a "step"

This module helps to compose flows and to construct the appropriate
shell/Hadoop/EMR scripts for running the flow. 
"""

import os
import hadoop

class Aggregation(object):
    """ Encapsulates information about the aggregation before a reducer. """
    
    def __init__(self, ntasks, ntasksPerReducer, nbin, nsort):
        self.ntasks = ntasks
        self.ntasksPerReducer = ntasksPerReducer
        assert nsort >= nbin
        self.nbin = nbin
        self.nsort = nsort
    
    def toLocalArgs(self, nReducers):
        args = []
        ntasks = self.ntasks
        if ntasks is None:
            ntasks = self.ntasksPerReducer * nReducers
        args.append('    --num-tasks=%d \\' % ntasks)
        args.append('    --bin-fields=%d \\' % self.nbin)
        args.append('    --sort-fields=%d \\' % self.nsort)
        return args
    
    def toHadoopArgs(self, config):
        begArgs, endArgs = [], []
        v = config.hadoopVersion
        confArg = hadoop.confArg(v)
        if self.ntasks is not None:
            ntasks = self.ntasks
        else:
            assert self.ntasksPerReducer is not None
            ntasks = self.ntasksPerReducer * config.nReducers
        begArgs.append('"%s", "mapred.reduce.tasks=%d",' % (confArg, ntasks))
        begArgs.append('"%s", "%s",' % (confArg, hadoop.keyFields(v, self.nbin)))
        endArgs.append('"-partitioner", "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",')
        begArgs.append('"%s", "stream.num.map.output.key.fields=%d",' % (confArg, self.nsort))
        return begArgs, endArgs

class PipelineConfig(object):
    """ Information about the particular system we're going to run on, which
        changes some of the exact arguments we use """
    
    def __init__(self, hadoopVersion, waitOnFail, emrStreamJar, nReducers, emrLocalDir):
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

class Step(object):
    """ Encapsulates a single step of the pipeline, i.e. a single
        MapReduce/Hadoop job """
    
    def __init__(\
        self,
        inps,
        output,
        name="(no name)",
        aggr=None,
        mapper=None,
        reducer=None,
        lineByLine=False,
        multipleOutput=False,
        cache=None,
        combiner=None):
        
        self.name = name
        self.inputs = inps
        self.output = output
        self.aggr = aggr
        self.mapper = mapper
        self.reducer = reducer
        self.lineByLine = lineByLine
        self.multipleOutput = multipleOutput
        self.cache = cache
        self.combiner = combiner
        self.inputFormat = None
        self.outputFormat = None
    
    def toLocalCmd(self, localConf):
        lines = []
        glue = None
        if self.mapper is not None:
            final = self.reducer is None
            lines.append('%PYPY% %BASE%/src/local/map.py \\')
            lines.append('    --name "%s" \\' % self.name)
            for inp in self.inputs:
                lines.append('    --input "%s" \\' % inp.toUrl())
            if final:
                if self.output is not None:
                    lines.append('    --output "%s" \\' % self.output.toUrl())
            else:
                # Pick an output directory that will feed reducer
                lines.append('    --output "%s" \\' % glue)
            if final and self.multipleOutput:
                lines.append('    --multiple-outputs \\')
            if self.lineByLine:   lines.append('    --line-by-line \\')
            if localConf.keepAll: lines.append('    --keep-all \\')
            if localConf.force:   lines.append('    --force \\')
            lines.append('    --num-processes %d \\' % localConf.numProcesses)
            # Separator between wrapper arguments and mapper command
            lines.append('    -- \\')
            lines.append('    ' + self.mapper)
            
            lines.append('''
if [ $? != 0 ] ; then
    echo "ERROR: map step %s failed with exitlevel $?"
    exit 1
fi
''' % (self.name))
        
        if self.reducer is not None:
            first = self.mapper is None
            lines.append('%PYPY% %BASE%/src/local/reduce.py \\')
            lines.append('    --name "%s" \\' % self.name)
            if first:
                lines.append('    --input "%s" \\' % ','.join([inp.toUrl() for inp in self.inputs]))
            else:
                lines.append('    --input "%s" \\' % glue)
            if self.output is not None:
                lines.append('    --output "%s" \\' % self.output.toUrl())
            if self.multipleOutput:
                lines.append('    --multiple-outputs \\')
            assert not first or not self.lineByLine
            if localConf.keepAll: lines.append('    --keep-all \\')
            if localConf.force:   lines.append('    --force \\')
            lines.append('    --num-processes %d \\' % localConf.numProcesses)
            lines.extend(self.aggr.toLocalArgs(localConf.numProcesses))
            # Separator between wrapper arguments and reducer command
            lines.append('    -- \\')
            lines.append('    ' + self.reducer)
            
            lines.append('''
if [ $? != 0 ] ; then
    echo "ERROR: reduce step %s failed with exitlevel $?"
    exit 1
fi
''' % (self.name))

        return '\n'.join(lines)
    
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
        
        begArgs, endArgs = [], []
        if self.aggr is not None:
            begArgs, endArgs = self.aggr.toHadoopArgs(config)
        
        if self.multipleOutput:
            self.outputFormat = 'edu.jhu.cs.MultipleOutputFormat'
            begArgs.append('"-libjars", "%BASE%/lib/multiplefiles.jar",')
        if self.lineByLine:
            self.inputFormat = 'org.apache.hadoop.mapred.lib.NLineInputFormat'
        
        endArgs.append('"-input", "%s",' % ','.join(map(lambda x: x.toUrl(), self.inputs)))
        if self.output is None:
            endArgs.append('"-output", "hdfs:///tmp/dummy1",')
        else:
            endArgs.append('"-output", "%s",' % self.output.toUrl())
        if self.mapper is not None:
            endArgs.append('"-mapper", "%s",' % self.mapper)
        else:
            endArgs.append('"-mapper", "cat",')

        if self.cache is not None:
            begArgs.append('"-archives", "%s",' % self.cache.toNativeUrl())
        
        if self.reducer is not None:
            endArgs.append('"-reducer", "%s",' % self.reducer)
        else:
            begArgs.append('"%s", "mapred.reduce.tasks=0",' % hadoop.confArg(config.hadoopVersion))
        
        if self.combiner is not None:
            endArgs.append('"-combiner", "%s",' % self.combiner)

        if self.inputFormat is not None:
            endArgs.append('"-inputformat", "%s",' % self.inputFormat)
        if self.outputFormat is not None:
            endArgs.append('"-outputformat", "%s",' % self.outputFormat)
        
        for a in begArgs + endArgs: lines.append('      ' + a)
        
        # Remove trailing comma
        lines[-1] = lines[-1][:-1]
        
        # Add cache files for any tools used by this job
        
        lines.append('    ]')
        lines.append('  }')
        lines.append('}')
        
        return '\n'.join(lines)

