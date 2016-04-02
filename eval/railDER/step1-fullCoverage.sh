#!/bin/sh


## Usage
# sh step1-fullCoverage.sh railGEU
# sh step1-fullCoverage.sh resub


# Define variables
EXPERIMENT=$1
SHORT="fullCov-${EXPERIMENT}"
CORES=10

# Directories
ROOTDIR=/dcl01/lieber/ajaffe/derRuns/railDER
MAINDIR=${ROOTDIR}/${EXPERIMENT}
WDIR=${MAINDIR}/CoverageInfo

if [[ "${EXPERIMENT}" == "railGEU" ]]
then
    DATADIR=/dcl01/lieber/ajaffe/derRuns/railDER/bigwig
    CUTOFF=5
elif [[ "${EXPERIMENT}" == "resub" ]]
then
    DATADIR=/dcl01/leek/data/geuvadis_rail_v0.1.9/coverage_bigwigs
    CUTOFF=5
else
    echo "Specify a valid experiment: railGEU, resub"
fi



# Construct shell file
echo 'Creating script for loading the Coverage data'
cat > ${ROOTDIR}/.${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=20G,h_vmem=40G,h_fsize=50G
#$ -N ${SHORT}
#$ -pe local ${CORES}

echo '**** Job starts ****'
date

# Make logs directory
mkdir -p ${WDIR}/logs

# Load the data, save the coverage without filtering, then save each file separately
cd ${WDIR}
module load R/3.1.x
Rscript ${ROOTDIR}/step1-fullCoverage.R -d "${DATADIR}" -p "bw" -c "${CUTOFF}" -m ${CORES}

## Move log files into the logs directory
mv ${ROOTDIR}/${SHORT}.* ${WDIR}/logs/

echo '**** Job ends ****'
date
EOF

call="qsub .${SHORT}.sh"
echo $call
$call
