"""
aws.py

Routines and variables for
"""

instNcores = {
    "m1.small"    : 1,
    "m1.large"    : 2,
    "m1.xlarge"   : 4,
    "c1.medium"   : 2,
    "c1.xlarge"   : 8,
    "m2.xlarge"   : 2,
    "m2.2xlarge"  : 4,
    "m2.4xlarge"  : 8,
    "cc1.4xlarge" : 8
}

def validateInstType(s): return s in instNcores

instSwap = {
    "m1.small"    : (2 *1024), #  1.7 GB
    "m1.large"    : (8 *1024), #  7.5 GB
    "m1.xlarge"   : (16*1024), # 15.0 GB
    "c1.medium"   : (2 *1024), #  1.7 GB
    "c1.xlarge"   : (8 *1024), #  7.0 GB
    "m2.xlarge"   : (16*1024), # 17.1 GB
    "m2.2xlarge"  : (16*1024), # 34.2 GB
    "m2.4xlarge"  : (16*1024), # 68.4 GB
    "cc1.4xlarge" : (16*1024)  # 23.0 GB
}

instBits = {
    "m1.small"    : 32,
    "m1.large"    : 64,
    "m1.xlarge"   : 64,
    "c1.medium"   : 32,
    "c1.xlarge"   : 64,
    "m2.xlarge"   : 64,
    "m2.2xlarge"  : 64,
    "m2.4xlarge"  : 64,
    "cc1.4xlarge" : 64
}

class EmrCluster(object):
    """ Encapsulates information about the nodes in the cluster """
    
    def __init__(self, itMaster, numMaster, itCore, numCore, itTask, numTask, bidPrice):
        assert numMaster > 0 and numCore > 0
        assert itMaster is not None and itMaster in instBits
        assert itCore is not None and itCore in instBits
        assert numTask == 0 or bidPrice is not None
        assert numTask == 0 or bidPrice > 0.0
        
        self.itMaster = itMaster
        self.itCore = itCore
        self.itTask = itTask
        self.numMaster = numMaster
        self.numCore = numCore
        self.numTask = numTask
        self.bidPrice = bidPrice
        
        # Get bit width on desired instances (32- or 64-)
        self.itBitsMaster, self.itBitsCore = instBits[itMaster], instBits[itCore]
        self.itBitsTask = 0
        if itTask is not None:
            self.itBitsTask = instBits[itTask]
        
        # Get amt of swap to add to desired instances
        self.itSwapMaster = instSwap[itMaster]
        self.itSwapCore = instSwap[itCore]
        self.itSwapTask = 0
        if itTask is not None:
            self.itSwapTask = instSwap[itTask]
    
    def numCores(self):
        return instNcores[self.itCore]
    
    def emrArgs(self):
        ret = []
        ret.append("--instance-group master --instance-type %s --instance-count %d" % (self.itMaster, self.numMaster))
        ret.append("--instance-group core --instance-type %s --instance-count %d" % (self.itCore, self.numCore))
        if self.numTask > 0:
            ret.append("--instance-group task --instance-type %s --instance-count %d --bid-price %f" % (self.itTask, self.numTask, self.bidPrice))
        return " ".join(ret)

def credsFromEnvironment():
    """ Return AWS credentials (AWS ID and AWS key) by parsing standard
        environment variables """
    import os
    g = os.environ.get
    return (g('AWS_ACCESS_KEY_ID') or os.environ.get('AWS_ACCESS_ID'),
            g('AWS_SECRET_ACCESS_KEY') or os.environ.get('AWS_ACCESS_KEY'))

def bootstrapHadoop(ncores=None, swap=None, memoryIntensive=False):
    """ Return Hadoop-related bootstrap actions for EMR cluster.  These can
        then be passed on the command line to the elastic-mapreduce script. """
    
    ret = []
    
    if memoryIntensive:
        # Memory-intensive mode
        ret.append('--bootstrap-action s3://elasticmapreduce/bootstrap-actions/configurations/latest/memory-intensive')
        ret.append('--bootstrap-name "Set memory-intensive mode"')
    
    # Don't reuse JVMs as much; this minimizes heap craziness
    ret.append('--bootstrap-action s3://elasticmapreduce/bootstrap-actions/configure-hadoop')
    ret.append('--bootstrap-name "Configure Hadoop"')
    #ret.append('--args "-s,mapred.job.reuse.jvm.num.tasks=1,-m,mapred.tasktracker.reduce.tasks.maximum=$cores,-s,io.sort.mb=100"')
    coresStr = ""
    if ncores is not None:
        coresStr = "-m,mapred.tasktracker.reduce.tasks.maximum=%n," % ncores
    ret.append('--args "%s-s,io.sort.mb=100"' % coresStr)
    
    if swap is not None:
        ret.append('--bootstrap-action s3://elasticmapreduce/bootstrap-actions/add-swap')
        ret.append('--bootstrap-name "Add Swap"')
        ret.append('--args "%s"' % swap)
    
    return ' '.join(ret)

def bootstrapFetch(name, url, emrLocalDir):
    """ Create a bootstrap action that copies the file at a given path within
        the given bucket into the given directory on the node. """
    ret = ['--bootstrap-action s3://tornado-emr/bootstrap/s3cmd_s3_tarball.sh',
           '--bootstrap-name "%s"' % name,
           '--args "%s,%s"' % (url.toNonNativeUrl(), emrLocalDir)]
    return ' '.join(ret)
