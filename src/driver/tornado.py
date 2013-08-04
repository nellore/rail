"""
tornado.py

Ben Langmead, 7/28/2013

Driver script for the Tornado pipeline.  Right now, just does Amazon Elastic
MapReduce mode.  Uses Amazon's elastic-mapreduce Ruby script to actually
launch the job.  This script creates a shell script and accompanying JSON file
that are used to run elastic-mapreduce.

  Preprocess
      |
    Align
  /             \
Merge          Intron
 |
Walk-prenorm
 |
Normalize
 |
Normalize-post

Walk-fit
 |
Ebayes
 |
HMM-params
 |
HMM
 |
aggr_path

TODO: Think about all the file plumbing:
- Tornado scripts
- Bowtie
- Reference jar (derived from iGenomes)

These all need to go somewhere.

"""

import os
import argparse
import string
import sys
import tempfile
import site

import aws
import pipeline
import tornado_pipeline
import tools
import tornado_config

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "util"))

import url
import path

parser = argparse.ArgumentParser(description='Generate and run a script for Tornado.')

#
# Modes of operation
#
parser.add_argument(\
    '--local', action='store_const', const=True, help='Run in local mode.')
parser.add_argument(\
    '--hadoop', action='store_const', const=True, help='Run in Hadoop mode.')
parser.add_argument(\
    '--emr', action='store_const', const=True, help='Run in Elastic MapReduce mode.')
parser.add_argument(\
    '--test', action='store_const', const=True, help='Just test to see if requisite files are around.')

#
# Basic plumbing
#
parser.add_argument(\
    '--manifest', metavar='PATH', type=str, required=False, help='URL for manifest file')
parser.add_argument(\
    '--input', metavar='PATH', type=str, required=False, help='URL for input directory')
parser.add_argument(\
    '--output', metavar='PATH', type=str, required=True, help='URL for output directory')
parser.add_argument(\
    '--reference', metavar='PATH', type=str, required=False, help='URL for reference archive')
parser.add_argument(\
    '--intermediate', metavar='PATH', type=str, help='URL for intermediate files')
parser.add_argument(\
    '--keep-intermediates', action='store_const', const=True, help='Keep intermediate files in addition to final output files.')
parser.add_argument(\
    '--dry-run', action='store_const', const=True, help='Just generate script for launching EMR cluster, but don\'t launch it.')

#
# Advanced plumbing
#
parser.add_argument(\
    '--just-preprocess', action='store_const', const=True, help='Just run the preprocessing step pipeline.')
parser.add_argument(\
    '--just-align', action='store_const', const=True, help='Just run the align step pipeline.')
parser.add_argument(\
    '--just-coverage', action='store_const', const=True, help='Just run the coverage step pipeline.')
parser.add_argument(\
    '--just-junction', action='store_const', const=True, help='Just run the junction step pipeline.')
parser.add_argument(\
    '--just-differential', action='store_const', const=True, help='Just run the differential pipeline.')
parser.add_argument(\
    '--start-with-align', action='store_const', const=True, help='Resume from just before the align pipeline.')
parser.add_argument(\
    '--start-with-coverage', action='store_const', const=True, help='Resume from just before the coverage pipeline.')
parser.add_argument(\
    '--start-with-junction', action='store_const', const=True, help='Resume from just before the junction pipeline.')
parser.add_argument(\
    '--start-with-differential', action='store_const', const=True, help='Resume from just before the differential pipeline.')
parser.add_argument(\
    '--stop-after-preprocess', action='store_const', const=True, help='Stop after the preprocessing pipeline.')
parser.add_argument(\
    '--stop-after-align', action='store_const', const=True, help='Stop after the align pipeline.')
parser.add_argument(\
    '--stop-after-coverage', action='store_const', const=True, help='Stop after the coverage pipeline.')
parser.add_argument(\
    '--stop-after-junction', action='store_const', const=True, help='Stop after the junction pipeline.')
parser.add_argument(\
    '--no-junction', action='store_const', const=True, help='Don\'t run the junction pipeline.')

#
# Hadoop params
#
parser.add_argument(\
    '--hadoop-script', metavar='PATH', type=str, help='Location of Hadoop script')
parser.add_argument(\
    '--streaming-jar', metavar='PATH', type=str, help='Location of Hadoop streaming jar')

#
# Elastic MapReduce params
#
parser.add_argument(\
    '--emr-script', metavar='PATH', type=str, help='Path to Amazon elastic-mapreduce script')
parser.add_argument(\
    '--emr-local-dir', metavar='PATH', type=str, help='Path to a local directory on the EMR nodes where the reference archive and Tornado scripts will be copied.')
parser.add_argument(\
    '--hadoop-version', metavar='VERS', type=str, help='Hadoop version number to use')
parser.add_argument(\
    '--credentials', metavar='PATH', type=str, help='Amazon Elastic MapReduce credentials file')
parser.add_argument(\
    '--alive', action='store_const', const=True, help='Keep cluster alive after it finishes job.')
parser.add_argument(\
    '--name', metavar='STR', type=str, help='Amazon Elastic MapReduce job name.')
parser.add_argument(\
    '--no-emr-debug', action='store_const', const=True, help='Don\'t enable EMR debugging functions.')
parser.add_argument(\
    '--instance-type', metavar='STR', type=str, help='Amazon EC2 instance type to use for all nodes.')
parser.add_argument(\
    '--instance-types', metavar='STR,STR,STR', type=str, help='Comma-separated list of EC2 instance types to use for MASTER, CORE and TASK nodes.')
parser.add_argument(\
    '--instance-count', metavar='INT', type=int, help='Number of EC2 instances to use (1 MASTER, rest will be CORE).')
parser.add_argument(\
    '--instance-counts', metavar='INT,INT,INT', type=str, help='Comma-separated list of number of EC2 instances to use for MASTER, CODE and TASK nodes.')
parser.add_argument(\
    '--bid-price', metavar='FLOAT', type=float, help='EC2 spot instance bid price.  If this is specified, TASK nodes will be spot instances.')

#
# Preprocessing params
#
parser.add_argument(\
    '--preprocess-output', metavar='PATH', type=str, help='Put output from preprocessing step here')
parser.add_argument(\
    '--preprocess-compress', metavar='gzip|bzip2', type=str, help='Type of compression to use for preprocessing output.')

tornado_config.addArgs(parser)

#
# Other params
#
parser.add_argument(\
    '--verbose', action='store_const', const=True, help='Print lots of info to stderr.')
parser.add_argument(\
    '--version', action='store_const', const=True, help='Just print version information and quit.')
parser.add_argument(\
    '--set-version', metavar='VER', type=str, help='Force Tornado to use a particular version.')

args = parser.parse_args()
if args.local:
    raise RuntimeError("--local mode not yet implemented")
if args.hadoop:
    raise RuntimeError("--hadoop mode not yet implemented")

def parseVersion():
    tornado_path = os.path.dirname(base_path)
    for basename in [ "VERSION", "TORNADO_VERSION" ]:
        vpath = os.path.join(tornado_path, basename)
        if os.path.exists(vpath):
            with open(vpath, 'r') as fh:
                return fh.readline().rstrip()
    else:
        raise RuntimeError("Could not find VERSION file (looked in '%s') and --set-version not specified" % tornado_path)

ver = args.set_version
if ver is None:
    ver = parseVersion()

print >> sys.stderr, "Tornado v" + ver

assert ver is not None

appName = "tornado" # app name

# TODO: also support local and Hadoop modes
mode = "emr"

# For URLs, remove trailing slash(es) and parse
inp = None
if args.input is not None:
    inp = url.Url(args.input.rstrip('/'))
manifest = None
if args.manifest is not None:
    manifest = url.Url(args.manifest.rstrip('/'))
out = url.Url(args.output.rstrip('/'))
logUrl = url.Url(out.toUrl() + '/' + "logs")
if args.intermediate is not None:
    intermediate = url.Url(args.intermediate.rstrip('/'))
else:
    intermediate = url.Url("hdfs://%s/intermediate" % appName)
ref = None
if args.reference is not None:
    ref = url.Url(args.reference.rstrip('/'))

jobName = None
swapAdd = 0
hadoopVerion, hadoopVersionToks = None, []
cred = None
emrScript = None
emrArgs = []
failAction = None
emrLocalDir = None
waitFail = False
emrCluster = None
emrStreamJar = None

#
# Parse AWS-related arguments.  Exceptions are raised when there are obvious
# issues.
#
if mode == 'emr':
    itMaster, itCore, itTask = None, None, None
    numMaster, numCore, numTask = 0, 0, 0
    
    # Parse instance type arguments
    if args.instance_types is not None:
        if args.instance_types.count(',') != 2:
            raise RuntimeError("Could not parse --instance-types string: '%s'" % args.instance_types)
        itMaster, itCore, itTask = string.split(args.instance_types, ',')
        if not aws.validateInstType(itMaster):
            raise RuntimeError("Invalid instance type for MASTER in --instance-types: '%s'" % itMaster)
        if not aws.validateInstType(itCore):
            raise RuntimeError("Invalid instance type for CODE in --instance-types: '%s'" % itCore)
        if not aws.validateInstType(itTask):
            raise RuntimeError("Invalid instance type for TASK in --instance-types: '%s'" % itTask)
    elif args.instance_type is not None:
        itMaster = args.instance_type
        itCore = args.instance_type
        itTask = None
        if not aws.validateInstType(itMaster):
            raise RuntimeError("Invalid instance type in --instance-type: '%s'" % itMaster)
    else:
        itMaster, itCore, itTask = "c1.xlarge", "c1.xlarge", None
    
    # Parse # instance arguments
    if args.instance_counts:
        c = args.instance_counts
        if c.count(',') != 2:
            raise RuntimeError("Could not parse --instance-counts string: '%s'" % c)
        numMaster, numCore, numTask = string.split(c, ',')
        numMaster, numCore, numTask = int(numMaster), int(numCore), int(numTask)
        usingSpot = numTask > 0
    elif args.instance_count:
        numMaster = 1
        numCore = int(args.instance_count)
    else:
        numMaster, numCore, numTask = 1, 1, 0
    
    emrCluster = aws.EmrCluster(itMaster, numMaster, itCore, numCore, itTask, numTask, args.bid_price)
    
    # Parse a few more EMR arguments
    jobName = args.name or "Tornado job"
    
    # Get AWS credentials
    cred = args.credentials
    if cred is not None and not path.is_exe(cred):
        raise RuntimeError("File specified with --credentials ('%s') doesn't exist or isn't executable" % args.credentials)
    
    # Get EMR script
    emrScript = args.emr_script
    if emrScript is not None and not path.is_exe(emrScript):
        raise RuntimeError("File specified with --emr-script ('%s') doesn't exist or isn't executable" % args.emrScript)
    
    # Sanity check Hadoop version
    # http://docs.aws.amazon.com/ElasticMapReduce/latest/DeveloperGuide/emr-plan-hadoop-version.html
    hadoopVersion = args.hadoop_version or "1.0.3"
    vers = map(int, string.split(hadoopVersion, '.'))
    if vers[0] < 1 and vers[1] < 20:
        raise RuntimeError("Not compatible with Hadoop versions before 0.20.  --hadoop-version was '%s'" % hadoopVersion)
    elif len(vers) < 2 or len(vers) > 4:
        raise RuntimeError("Could not parse --hadoop-version '%s'" % hadoopVersion)
    elif len(vers) >= 3 and vers[:3] == [1, 0, 3]:
        emrArgs.append("--hadoop-version=1.0.3")
        emrArgs.append("--ami-version 2.4")
    elif len(vers) >= 3 and vers[:3] == [0, 20, 205]:
        emrArgs.append("--hadoop-version=0.20.205")
        emrArgs.append("--ami-version 2.0")
    elif vers[:2] == [0, 20]:
        emrArgs.append("--hadoop-version=0.20")
        emrArgs.append("--ami-version 1.0")
    else:
        raise RuntimeError("Unexpected --hadoop-version '%s' (%s)" % (hadoopVersion, str(vers)))
    hadoopVersionToks = vers
    
    # Suppress debugging?
    if not args.no_emr_debug:
        emrArgs.append("--enable-debugging")
    
    # ActionOnFailure
    waitFail = args.alive
    
    # Local dir, where ref genome and Tornado scripts are installed
    emrLocalDir = args.emr_local_dir or "/mnt"
    
    # Set up name of streaming jar for EMR mode
    emrStreamJar = "/home/hadoop/contrib/streaming/hadoop-streaming-%s.jar" % hadoopVersion
    if hadoopVersionToks[0] == 0 and hadoopVersionToks[1] < 205:
        emrStreamJar = "/home/hadoop/contrib/streaming/hadoop-%s-streaming.jar" % hadoopVersion
    
    # Sanity-check URLs
    if inp is not None and not inp.isS3():
        raise RuntimeError("--input argument '%s' is not an S3 URL" % inp)
    if not out.isS3():
        raise RuntimeError("--output argument '%s' is not an S3 URL" % out)
    if ref is not None and not ref.isS3():
        raise RuntimeError("--reference argument '%s' is not an S3 URL" % ref)
    if args.intermediate is not None and not intermediate.isS3():
        raise RuntimeError("--intermediate argument '%s' is not an S3 URL" % args.intermediate)

tconf = tornado_config.TornadoConfig(args)
pconf = pipeline.PipelineConfig(hadoopVersionToks, waitFail, emrStreamJar, emrCluster.numCoreProcessors(), emrLocalDir, args.preprocess_compress, out)

pipelines = ["preprocess", "align", "coverage", "junction", "differential"]

# Might start partway through
if args.start_with_align:
    pipelines.remove("preprocess")
if args.start_with_junction or args.start_with_coverage:
    for s in ["preprocess", "align"]: pipelines.remove(s)
if args.start_with_differential:
    for s in ["preprocess", "align", "junction", "coverage"]: pipelines.remove(s)

# Might end partway through
if args.stop_after_preprocess:
    for s in ["align", "junction", "coverage", "differential"]: pipelines.remove(s)
if args.stop_after_align:
    for s in ["junction", "coverage", "differential"]: pipelines.remove(s)
if args.stop_after_coverage or args.stop_after_junction:
    pipelines.remove("differential")

# Might skip junction pipeline
if args.no_junction:
    pipelines.remove("junction")

# Might just be running one pipeline
if args.just_preprocess:
    pipelines = ["preprocess"]
elif args.just_align:
    pipelines = ["align"]
elif args.just_junction:
    pipelines = ["junction"]
elif args.just_coverage:
    pipelines = ["coverage"]
elif args.just_differential:
    pipelines = ["differential"]
assert len(pipelines) > 0

pipelineSteps = {
    'preprocess'   : ['preprocess'],
    'align'        : ['align'],
    'junction'     : ['intron'],
    'coverage'     : ['merge', 'walk_prenorm', 'normalize', 'normalize_post'],
    'differential' : ['walk_fit', 'ebayes', 'hmm_params', 'hmm', 'aggr_path'] }

useBowtie = 'align' in pipelines
useIndex = 'align' in pipelines
useGtf = 'align' in pipelines
useFasta = 'align' in pipelines or 'junction' in pipelines
useRef = useIndex or useGtf or useFasta
useKenttools = 'coverage' in pipelines or 'differential' in pipelines
useSamtools = False
useSraToolkit = 'preprocess' in pipelines
useR = 'differential' in pipelines
useManifest = 'coverage' in pipelines
useInput = 'preprocess' not in pipelines

if useManifest and manifest is None:
    raise RuntimeError("Must specify --manifest when job involves preprocessing")

if useRef and ref is None:
    raise RuntimeError("Must specify --reference when job involves alignment or junction pipelines")

if useInput and inp is None:
    raise RuntimeError("Must specify --input when job does not involve preprocessing")

allSteps = [ None ] + [ i for sub in map(pipelineSteps.get, pipelines) for i in sub ] + [ None ]

stepClasses = {\
    'preprocess'     : tornado_pipeline.PreprocessingStep,
    'align'          : tornado_pipeline.AlignStep,
    'intron'         : tornado_pipeline.IntronStep,
    'merge'          : tornado_pipeline.MergeStep,
    'walk_prenorm'   : tornado_pipeline.WalkPrenormStep,
    'normalize'      : tornado_pipeline.NormalizeStep,
    'normalize_post' : tornado_pipeline.NormalizePostStep,
    'walk_fit'       : tornado_pipeline.WalkFitStep,
    'ebayes'         : tornado_pipeline.EbayesStep,
    'hmm_params'     : tornado_pipeline.HmmParamsStep,
    'hmm'            : tornado_pipeline.HmmStep,
    'aggr_path'      : tornado_pipeline.AggrPathStep }

inDirs, outDirs, steps = [], [], []
for prv, cur, nxt in [ allSteps[i:i+3] for i in xrange(0, len(allSteps)-2) ]:
    assert cur in stepClasses
    if prv is None:
        if cur == "preprocess":
            inDirs.append(manifest.toUrl())
        else:
            inDirs.append(inp.toUrl())
    else: inDirs.append(outDirs[-1])
    if nxt is None: outDirs.append(out.toUrl() + "/final")
    else:
        if args.preprocess_output and cur == "preprocess":
            outDirs.append(args.preprocess_output)
        else:
            outDirs.append("%s/%s_out" % (intermediate.toUrl(), cur))
    steps.append(stepClasses[cur](inDirs[-1], outDirs[-1], tconf, pconf))

if mode == 'emr':
    jsonStr = "[\n" + "\n".join([ step.toEmrCmd(pconf) for step in steps ]) + "\n]\n"
    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as jsonFh:
        jsonFn = jsonFh.name
        jsonFh.write(jsonStr)
    
    print >> sys.stderr, "Json output in: '%s'" % jsonFn
    
    if emrScript is None:
        emrScript = "elastic-mapreduce"
    
    cmdl = [
        emrScript, "--create", # create job flow
        "--name", jobName,     # name job flow
        "--json", jsonFn,
        "--log-uri", logUrl.toUrl() ]     # file with JSON description of flow

    if cred is not None:
        cmdl.extend(["-c", cred])
    
    cmdl.append(emrCluster.emrArgs())
    tornadoUrl = url.Url("s3://tornado-emr/bin/tornado-%s.tar.gz" % ver)
    if useRef:
        cmdl.append(aws.bootstrapFetchTarball("Fetch Tornado ref archive", ref, emrLocalDir))
    if useManifest:
        cmdl.append(aws.bootstrapFetchFile("Fetch manifest file", manifest, emrLocalDir, "MANIFEST"))
    cmdl.append(tools.bootstrapTool("python"))
    if useBowtie:
        cmdl.append(tools.bootstrapTool("bowtie"))
    if useSraToolkit:
        cmdl.append(tools.bootstrapTool("sra-toolkit"))
    if useR:
        cmdl.append(tools.bootstrapTool("R"))
    if useKenttools:
        cmdl.append(tools.bootstrapTool("kenttools", dest="/mnt/bin"))
    if useSamtools:
        cmdl.append(tools.bootstrapTool("samtools"))
    # Get Tornado scripts and run Makefile for swig code
    cmdl.append(tools.bootstrapTool("tornado", src=tornadoUrl, dest="/mnt"))
    cmdl.extend(emrArgs)
    
    cmd = ' '.join(cmdl)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".sh", delete=False) as shFh:
        shFn = shFh.name
        shFh.write(cmd + '\n')
    print >> sys.stderr, "elastic-mapreduce command in: '%s'" % shFn
    if not args.dry_run:
        os.system(cmd)
    else:
        print >> sys.stderr, "--dry-run mode: not running elastic-mapreduce script"
