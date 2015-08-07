#!/bin/sh

## Usage
# sh step6b-regionMatrix-merge.sh railGEU
# sh step6b-regionMatrix-merge.sh resub

# Define variables
EXPERIMENT=$1
SHORT="regMat-merge-${EXPERIMENT}"

# Directories
ROOTDIR=/dcs01/ajaffe/Brain/derRuns/railDER
MAINDIR=${ROOTDIR}/${EXPERIMENT}
WDIR=${MAINDIR}/regionMatrix

if [[ "${EXPERIMENT}" == "railGEU" ]]
then
    CUTOFF=5
elif [[ "${EXPERIMENT}" == "resub" ]]
then
    CUTOFF=5
else
    echo "Specify a valid experiment: railGEU, resub"
fi

# Construct shell file
sname="${SHORT}"
echo 'Creating script for loading the Coverage data'
cat > ${ROOTDIR}/.${SHORT}.sh <<EOF
#!/bin/bash	
#$ -cwd
#$ -m e
#$ -l mem_free=75G,h_vmem=180G,h_fsize=40G
#$ -N ${sname}
#$ -hold_jid regMat-${EXPERIMENT}-chr*

echo "**** Job starts ****"
date

# Make logs directory
mkdir -p ${WDIR}/logs

# Load coverage & get region matrix
cd ${WDIR}
module load R/3.1.x
R -e "cutoff <- ${CUTOFF}; source('${ROOTDIR}/step6b-regionMatrix-merge.R')"

## Move log files into the logs directory
mv ${ROOTDIR}/${sname}.* ${WDIR}/logs/

echo "**** Job ends ****"
date
EOF

call="qsub .${SHORT}.sh"
echo $call
$call
