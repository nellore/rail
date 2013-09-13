"""
flux.py

Wrapper for flux simulator that generates the .par file automatically from
various command-line options.  Also potentially handles many technical
replicates in parallel.
"""

import argparse
import subprocess
from multiprocessing import Pool

parser = argparse.ArgumentParser(description=\
                                     'Creates simulated data via flux simulator.')
parser.add_argument(\
    '--name', type=str, required=True, help='Name of the experiment.  All the flux files will have this as a prefix.'
)
parser.add_argument(\
    '--output-dir',type=str,required=True,help='Output directory of simulated reads'
)
parser.add_argument(\
    '--flux-path',type=str,required=True,help='Path of flux-simulator home directory'
)
parser.add_argument(\
    '--gtf-file',type=str,required=True,help='The fixed gtf file'
)
parser.add_argument(\
    '--chromosomes',type=str,required=True,help='The chromosomes directory'
)
parser.add_argument(\
    '--num-processes', type=int, default=1, help='Number of flux-simulator processes are allowed to run at once'
)
parser.add_argument(\
    '--num-samples',type=int,required=True,help='Number of simulated samples'
)
parser.add_argument(\
    '--num-reads',type=int,required=False,default=500000,help='Number of simulated reads'
)
parser.add_argument(\
    '--num-molecules',type=int,required=False,default=500,help='Number of simulated molecules'
)
parser.add_argument(\
    '--tss-mean',type=int,required=False,default=50,help='Average deviation from the annotated transcription start site (TSS)'
)
parser.add_argument(\
    '--polya-scale',type=str,required=False,default='NaN',help='Scale of the Weibull distribution, shifts the average length of poly-A tail sizes'
)
parser.add_argument(\
    '--polya-shape',type=str,required=False,default='NaN',help='Shape of the Weibull distribution describing poly-A tail sizes'
)
parser.add_argument(\
    '--frag-ur-eta',type=int,required=False,default=350,help='Average expected framgent size after fragmentations, i.e., number of breaks per unit length (exhautiveness of fragmentation)'
)
parser.add_argument(\
    '--frag-ur-d0',type=int,required=False,default=1,help='Minimum length of fragments produced by UR fragmentation'
)
parser.add_argument(\
    '--rt-min',type=int,required=False,default=500,help='Minimum length observed after reverse transcription of full-length transcripts'
)
parser.add_argument(\
    '--rt-max',type=int,required=False,default=5500,help='Maximum length observed after reverse transcription of full-length transcripts'
)
parser.add_argument(\
    '--pcr-probability',type=float,required=False,default=0.05,help='Default PCR distribution with 15 rounds and 20 bins'
)
parser.add_argument(\
    '--read-length',type=int,required=False,default=76,help='Length of reads'
)
parser.add_argument(\
    '--err-file',type=int,required=False,default=76,help='Error file specifications'
)
args=parser.parse_args()

def createParameterFile(par_name):
    """ Creates parameter file for flux-simulator.  See:
        http://sammeth.net/confluence/display/SIM/.PAR+Simulation+Parameters
    """
    par_out = open("%s/%s"%(args.output_dir,par_name),'w')
    par_out.write("NB_MOLECULES\t%d\n"%args.num_molecules)
    par_out.write("REF_FILE_NAME\t%s\n"%args.gtf_file)
    par_out.write("GEN_DIR\t%s\n"%args.chromosomes)
    par_out.write("LOAD_NONCODING\tNO\n")
    par_out.write("TSS_MEAN\t%d\n"%args.tss_mean)
    par_out.write("POLYA_SCALE\t%s\n"%args.polya_scale)
    par_out.write("POLYA_SHAPE\t%s\n"%args.polya_shape)
    par_out.write("FRAG_SUBSTRATE\tRNA\n")
    par_out.write("FRAG_METHOD\tUR\n")
    par_out.write("FRAG_UR_ETA\t%d\n"%args.frag_ur_eta)
    par_out.write("FRAG_UR_D0\t%d\n"%args.frag_ur_d0)
    par_out.write("RTRANSCRIPTION\tYES\n")
    par_out.write("RT_PRIMER\tRH\n")
    par_out.write("RT_LOSSLESS\tYES\n")
    par_out.write("RT_MIN\t%d\n"%args.rt_min)
    par_out.write("RT_MAX\t%d\n"%args.rt_max)
    par_out.write("GC_MEAN\tNaN\n")
    par_out.write("PCR_PROBABILITY\t%f\n"% args.pcr_probability)
    par_out.write("FILTERING\tNO\n")
    par_out.write("READ_NUMBER\t%d\n"%args.num_reads)
    par_out.write("READ_LENGTH\t%d\n"%args.read_length)
    par_out.write("PAIRED_END\tYES\n")
    par_out.write("ERR_FILE\t%d\n"%args.err_file)
    par_out.write("FASTA\tYES\n")
    par_out.write("UNIQUE_IDS\tYES\n")
    par_out.close()

def runFlux(par_name):
    """ Runs Flux Simulator on a particular sample name """
    flux_cmd="%s/bin/flux-simulator -p %s/%s"%(args.flux_path,args.output_dir,par_name)
    flux_proc = subprocess.Popen(flux_cmd, bufsize=-1, shell=True)
    return flux_proc #return process for parallelism

def createManifest(manifest_name,samples):
    """ Create a flux-simulator parameters file """
    manifest_out = open("%s/%s"%(args.output_dir,manifest_name),'w')
    #TODO: Need to correct the manifest format below
    for samp in samples:
        manifest_out.write("%s\t0\t%s\n"%(samp,samp))
    manifest_out.close()

def doTechRep(nm):
    """ Run flux-simulator for a single technical replicate.  Argument is the
        sample name. """
    createParameterFile(nm + ".par")
    exitlevel = runFlux(nm + ".par").wait()
    if exitlevel != 0:
        raise RuntimeError("Got exitlevel %d from a flux-simulator process" % exitlevel)

if __name__=="__main__":
    pool = Pool(args.num_processes)
    pool.map(doTechRep, [ "%s-%d" % (args.name, i) for i in xrange(args.num_samples) ])
