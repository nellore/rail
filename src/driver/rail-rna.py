"""
rail-rna.py

Ben Langmead, 7/28/2013

Driver script for the Rail-RNA pipeline. Uses Amazon's elastic-mapreduce Ruby script
to actuallylaunch the job.  This script creates a shell script and accompanying
JSON file that are used to run elastic-mapreduce.

WRITE NEW DIAGRAM OF RAIL HERE


TODO: Think about all the file plumbing:
- Rail-RNA scripts
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
# hyphenated filename bars conventional import usages
rail_rna_pipeline = __import__('rail-rna_pipeline')
rail_rna_config = __import__('rail-rna_config')
from config import addConfigArgs, GenericConfig
from tools import ToolConfigLocal, ToolConfigHadoop, ToolConfigEmr, bootstrapTool
from ref import RefConfigLocal, RefConfigHadoop, RefConfigEmr
from app import AppConfigLocal, AppConfigHadoop, AppConfigEmr
from files import FileConfigLocal, FileConfigHadoop, FileConfigEmr, addFileArgs
from tempfile import mkdtemp
from emr_mode import addEmrModeArgs
from hadoop_mode import addHadoopModeArgs
from local_mode import addLocalModeArgs, LocalConfig

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, "util"))

from url import Url
from path import is_exe

parser = argparse.ArgumentParser(description='Generate and run a script for Rail-RNA.')

addConfigArgs(parser)
addFileArgs(parser)

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
    '--start-with-preprocess', action='store_const', const=True, help='Start pipeline from preprocessing step.')
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
    '--no-coverage', action='store_const', const=True, help='Don\'t run the coverage pipeline.')
parser.add_argument(\
    '--no-junction', action='store_const', const=True, help='Don\'t run the junction pipeline.')
parser.add_argument(\
    '--no-differential', action='store_const', const=True, help='Don\'t run the differential expression pipeline.')

addEmrModeArgs(parser)
addHadoopModeArgs(parser)
addLocalModeArgs(parser)

rail_rna_config.addArgs(parser)

args = parser.parse_args()
if args.hadoop:
    raise RuntimeError("--hadoop mode not yet implemented")

appName = "Rail-RNA" # app name

def parseVersion():
    path = os.path.dirname(base_path)
    for basename in [ "VERSION", appName.upper() + "_VERSION" ]:
        vpath = os.path.join(path, basename)
        if os.path.exists(vpath):
            with open(vpath, 'r') as fh:
                return fh.readline().rstrip()
    else:
        raise RuntimeError("Could not find VERSION file (looked in '%s') and --set-version not specified" % path)

ver = args.set_version
if ver is None:
    ver = parseVersion()

print >> sys.stderr, "Rail-RNA v" + ver

assert ver is not None

mode = "emr"
if args.local: mode = 'local'
if args.hadoop: mode = 'hadoop'

# For URLs, remove trailing slash(es) and parse
inp = None
if args.input is not None:
    inp = Url(args.input.rstrip('/'))
manifest = None
if args.manifest is not None:
    manifest = Url(args.manifest.rstrip('/'))
out = Url(args.output.rstrip('/'))
logUrl = Url(out.toUrl() + '/' + "logs")
intermediate = None
if args.intermediate is not None:
    intermediate = Url(args.intermediate.rstrip('/'))
reference = None
if args.reference is not None:
    reference = Url(args.reference.rstrip('/'))

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
        itTask = args.instance_type
        if not aws.validateInstType(itMaster):
            raise RuntimeError("Invalid instance type in --instance-type: '%s'" % itMaster)
    else:
        itMaster, itCore, itTask = "c1.xlarge", "c1.xlarge", "c1.xlarge"
    
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
    jobName = args.name or appName + " job"
    
    # Get AWS credentials
    cred = args.credentials
    if cred is not None and not is_exe(cred):
        raise RuntimeError("File specified with --credentials ('%s') doesn't exist or isn't executable" % args.credentials)
    
    # Get EMR script
    emrScript = args.emr_script
    if emrScript is not None and not is_exe(emrScript):
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
        emrArgs.append("--ami-version 2.4.2")
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
    if waitFail:
        emrArgs.append("--alive --with-termination-protection")
    
    # Local dir, where ref genome and scripts are installed
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
    if reference is not None and not reference.isS3():
        raise RuntimeError("--reference argument '%s' is not an S3 URL" % reference)
    if args.intermediate is not None and not intermediate.isS3():
        raise RuntimeError("--intermediate argument '%s' is not an S3 URL" % args.intermediate)

# Parameters governing the algorithm
tconf = rail_rna_config.Rail_RNAConfig(args)
tconf.keep_alive = (False if mode == 'local' else True)
gconf = GenericConfig(args, out, intermediate)

pipelines = set(["preprocess", "align", "align_out", "coverage", "differential"])

# Might start partway through
if args.start_with_align:
    pipelines = pipelines - set(["preprocess"])
if args.start_with_junction or args.start_with_coverage:
    pipelines = pipelines - set(["preprocess", "align"])
if args.start_with_differential:
    pipelines = pipelines - set(["preprocess", "align", "align_out", "coverage"])

# Might end partway through
if args.stop_after_preprocess:
    pipelines = pipelines - set(["align", "junction", "coverage", "differential"])
if args.stop_after_align:
    pipelines = pipelines - set(["junction", "coverage", "differential"])
if args.stop_after_coverage or args.stop_after_junction:
    pipelines = pipelines - set(["junction", "coverage", "differential"])

if args.no_coverage:
    pipelines = pipelines - set(["coverage"])
if args.no_junction:
    pipelines = pipelines - set(["junction"])
if args.no_differential:
    pipelines = pipelines - set(["differential"])

# Might just be running one pipeline
if args.just_preprocess:
    pipelines = set(["preprocess"])
elif args.just_align:
    pipelines = set(["align"])
elif args.just_coverage:
    pipelines = set(["coverage"])
elif args.just_differential:
    pipelines = set(["differential"])
assert len(pipelines) > 0

useBowtie = 'align' in pipelines
useIndex = 'align' in pipelines
useGtf = False and 'align' in pipelines
useFasta = 'align' in pipelines or 'junction' in pipelines
useRef = useIndex or useGtf or useFasta
useKenttools = 'coverage' in pipelines or 'differential' in pipelines
useSamtools = 'align' in pipelines
useSraToolkit = 'preprocess' in pipelines
useR = 'differential' in pipelines
useManifest = 'coverage' in pipelines
useInput = 'preprocess' not in pipelines

if useManifest and manifest is None:
    raise RuntimeError("Must specify --manifest when job involves the coverage pipeline")

if useRef and reference is None:
    raise RuntimeError("Must specify --reference when job involves alignment or junction pipelines")

if useInput and inp is None:
    raise RuntimeError("Must specify --input when job does not involve preprocessing")

pipelineSteps = {
    'preprocess'   : ['preprocess'],
    'align'        : ['align', 'intron', 'intron_config', 'intron_fasta', 'intron_index', 'realign'],
    'align_out'    : ['bam', 'bed_pre', 'bed'],
    #'coverage'     : ['normalize_pre', 'normalize'],#, 'normalize_post'],
    #'coverage'      : [],
    'coverage'     : ['collapse', 'coverage_pre', 'coverage', 'coverage_post'],
    'differential' : ['walk_fit', 'ebayes', 'hmm_params', 'hmm', 'aggr_path'] }

allSteps = [ i for sub in map(pipelineSteps.get, pipelines) for i in sub ]

stepInfo = {\
    'preprocess'     : ([                                       ], rail_rna_pipeline.PreprocessingStep),
    'align'          : ([('preprocess',     '/push'            )], rail_rna_pipeline.AlignStep),
    'intron'         : ([('align',          '/intron'          )], rail_rna_pipeline.IntronStep),
    'intron_post'    : ([('align',          '/max_len'         ),
                         ('intron',         '/intron'          )], rail_rna_pipeline.IntronPostStep),
    'intron_config'  : ([('align',          '/max_len'         ),
                         ('intron',         '/intron'          )], rail_rna_pipeline.IntronConfigStep),
    'intron_fasta'   : ([('intron_config',  '/intron'          )], rail_rna_pipeline.IntronFastaStep),
    'intron_index'   : ([('intron_fasta',   '/intron'          )], rail_rna_pipeline.IntronIndexStep),
    'bam'            : ([('align',          '/end_to_end_sam'  ),
                         ('realign',        '/splice_sam'      )], rail_rna_pipeline.BamStep),
    'bed_pre'        : ([('realign',        '/intron'          )], rail_rna_pipeline.BedPreStep),
    'bed'            : ([('bed_pre',        '/bed'             )], rail_rna_pipeline.BedStep),
    'realign'        : ([('align',          '/unmapped'        )], rail_rna_pipeline.RealignStep),
    'coverage_pre'   : ([('collapse',       '/collapsed'       )], rail_rna_pipeline.CoveragePreStep),
    'normalize_pre'  : ([('align',          '/exon_diff'       )], rail_rna_pipeline.NormalizePreStep),
    'normalize'      : ([('normalize_pre',  '/o'               )], rail_rna_pipeline.NormalizeStep),
    'collapse'       : ([('align',          '/exon_diff'       ), 
                         ('realign',        '/exon_diff'       )], rail_rna_pipeline.CollapseStep),
    'coverage'       : ([('coverage_pre',   '/coverage'        )], rail_rna_pipeline.CoverageStep),
    'coverage_post'  : ([('coverage',       ''                 )], rail_rna_pipeline.CoveragePostStep),
    'normalize_post' : ([('normalize',      ''                 )], rail_rna_pipeline.NormalizePostStep),
    'walk_fit'       : ([('normalize_post', ''                 )], rail_rna_pipeline.WalkFitStep),
    'ebayes'         : ([('walk_fit',       ''                 )], rail_rna_pipeline.EbayesStep),
    'hmm_params'     : ([('ebayes',         ''                 )], rail_rna_pipeline.HmmParamsStep),
    'hmm'            : ([('hmm_params',     ''                 )], rail_rna_pipeline.HmmStep),
    'aggr_path'      : ([('hmm',            ''                 )], rail_rna_pipeline.AggrPathStep) }

# 'normalize_post' sends pushes normalization-factor .tsv to out/normalization_factors.tsv
# 'normalize' sends per-sample coverage bigBed to out
# 'intron' sends junction into to out

def buildFlow(allSteps, stepInfo, inp, out, inter, manifest):
    steps = []
    for cur in allSteps:
        prvs, cl = stepInfo[cur]
        indirs = []
        if len(prvs) == 0: indirs.append(manifest)
        for prv, subdir in prvs:
            if prv in allSteps:
                indirs.append(inter.plus(prv + subdir))
            else:
                indirs.append(inp)
        outdir = inter.plus(cur)
        if cl is not rail_rna_pipeline.RealignStep:
            steps.append(cl(indirs, outdir, tconf, gconf))
        else:
            '''Must cache intron index: see 
            http://docs.aws.amazon.com/ElasticMapReduce/latest/
            DeveloperGuide/emr-plan-input-distributed-cache.html'''
            steps.append(cl(indirs, outdir, tconf, gconf, cache=out.plus('index/intron_index.tar.gz#intron_index')))
    return steps

if intermediate is None:
    if mode == 'local':
        intermediate = Url(os.path.join(mkdtemp(), 'intermediate'))
        if args.verbose:
            print >> sys.stderr, 'Temporary directory for local mode: "%s"' % intermediate
    elif mode == 'hadoop':
        intermediate = Url("hdfs:///%s/intermediate" % appName)
    elif mode == 'emr':
        intermediate = Url("hdfs:///%s/intermediate" % appName)
assert intermediate is not None

# Need to set default intermediate dir
steps = buildFlow(allSteps, stepInfo, inp, out, intermediate, manifest)

if mode == 'local':
    appConf = AppConfigLocal(appName, True)
    toolConf = ToolConfigLocal(appName, True)
    refConf, fileConf = None, None
    if useRef:
        refConf = RefConfigLocal(reference.toUrl(), intermediate.toUrl(), out.toUrl(), args.igenomes, True)
    if useManifest:
        fileConf = FileConfigLocal(manifest, True)
    localConf = LocalConfig(args)
    
    shStr = '#!/bin/sh\n'
    if localConf.force:
        shStr += 'rm -rf %s\n' % out.toUrl()
    shStr += '\n'.join([ step.toLocalCmd(localConf) for step in steps ])
    if refConf is not None: shStr = refConf.config(shStr)
    if fileConf is not None: shStr = fileConf.config(shStr)
    shStr = toolConf.config(appConf.config(shStr))
    
    with tempfile.NamedTemporaryFile(mode="w", suffix=".sh", delete=False) as shFh:
        shFn = shFh.name
        shFh.write(shStr)
    print >> sys.stderr, "Shell output in: '%s'" % shFn
    
    cmdl = [ 'sh', shFn ]
    cmd = ' '.join(cmdl)
    if not args.dry_run:
        os.system(cmd)
    else:
        print >> sys.stderr, "--dry-run mode: not running shell script"
elif mode == 'hadoop':
    raise RuntimeError('--hadoop mode not implemented yet')
elif mode == 'emr':
    
    appConf = AppConfigEmr(emrLocalDir)
    toolConf = ToolConfigEmr(emrLocalDir)
    refConf, fileConf = None, None
    if useRef: refConf = RefConfigEmr(emrLocalDir, args.igenomes)
    if useManifest: fileConf = FileConfigEmr(emrLocalDir)
    
    # Parameters governing how to run the pipeline
    pconf = pipeline.PipelineConfig(hadoopVersionToks, waitFail, emrStreamJar, emrCluster.numCoreOrTaskProcessors(), emrLocalDir)

    jsonStr = "[\n" + ",\n".join([ step.toEmrCmd(pconf) for step in steps ]) + "\n]\n"
    if refConf is not None: jsonStr = refConf.config(jsonStr)
    if fileConf is not None: jsonStr = fileConf.config(jsonStr)
    jsonStr = toolConf.config(appConf.config(jsonStr))
    
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
    rail_RNAUrl = Url("s3://tornado-emr/bin/rail-rna-%s.tar.gz" % ver)
    
    cmdl.append(bootstrapTool("python"))
    if useBowtie:
        cmdl.append(bootstrapTool("bowtie"))
    if useSraToolkit:
        cmdl.append(bootstrapTool("sra-toolkit"))
    if useR:
        cmdl.append(bootstrapTool("R"))
    if useKenttools:
        cmdl.append(bootstrapTool("kenttools", dest="/mnt/bin"))
    if useSamtools:
        cmdl.append(bootstrapTool("samtools"))
    # Get Rail-RNA scripts and run Makefile for swig code
    cmdl.append(bootstrapTool("rail", src=rail_RNAUrl, dest="/mnt"))
    tarballs = []
    if useIndex:
        tarballs.append(Url(reference.toUrl().replace('.tar.gz', '.index.tar.gz')))
    if useGtf:
        tarballs.append(Url(reference.toUrl().replace('.tar.gz', '.gtf.tar.gz')))
    if useFasta:
        tarballs.append(Url(reference.toUrl().replace('.tar.gz', '.fasta.tar.gz')))
    if len(tarballs) > 0:
        cmdl.append(aws.bootstrapFetchTarballs("reference archives", tarballs, emrLocalDir))
    if useManifest:
        cmdl.append(aws.bootstrapFetchFile("manifest file", manifest, emrLocalDir, "MANIFEST"))
    #if not args.no_memory_intensive:
    #    # Memory-intensive mode
    #    cmdl.append('--bootstrap-action')
    #    cmdl.append('s3://elasticmapreduce/bootstrap-actions/configurations/latest/memory-intensive')
    #    cmdl.append('--bootstrap-name')
    #    cmdl.append('"set memory-intensive"')
    if args.ganglia:
        cmdl.append('--bootstrap-action')
        cmdl.append('s3://elasticmapreduce/bootstrap-actions/install-ganglia')
        cmdl.append('--bootstrap-name')
        cmdl.append('"ganglia"')
    if not args.no_add_swap:
        cmdl.append('--bootstrap-action')
        cmdl.append('s3://elasticmapreduce/bootstrap-actions/add-swap')
        cmdl.append('--bootstrap-name')
        cmdl.append('"add swap"')
        cmdl.append('--args "%d"' % emrCluster.swap())
    numProcessorsPerCore = emrCluster.numProcessorsPerCoreInstance()
    hadoopConfigs = []
    hadoopConfigs.append('-s,mapred.job.reuse.jvm.num.tasks=1')
    hadoopConfigs.append('-s,mapred.tasktracker.reduce.tasks.maximum=%d' % numProcessorsPerCore)
    hadoopConfigs.append('-s,mapred.tasktracker.map.tasks.maximum=%d' % numProcessorsPerCore)
    cmdl.append('--bootstrap-action')
    cmdl.append('s3://elasticmapreduce/bootstrap-actions/configure-hadoop')
    cmdl.append('--bootstrap-name')
    cmdl.append('"configure hadoop"')
    if not                      args.enable_speculative:
        hadoopConfigs.append('-m,mapred.map.tasks.speculative.execution=false')
        hadoopConfigs.append('-m,mapred.reduce.tasks.speculative.execution=false')
    cmdl.append('--args "%s"' % ','.join(hadoopConfigs))
    
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
