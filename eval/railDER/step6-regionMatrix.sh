#!/bin/sh

## Usage
# sh step6-regionMatrix.sh railGEU
# sh step6-regionMatrix.sh resub

# Define variables
EXPERIMENT=$1
SHORT="regMat-${EXPERIMENT}"

# Directories
ROOTDIR=/dcl01/lieber/ajaffe/derRuns/railDER
MAINDIR=${ROOTDIR}/${EXPERIMENT}
WDIR=${MAINDIR}/regionMatrix

if [[ "${EXPERIMENT}" == "railGEU" ]]
then
    CUTOFF=5
    RLENGTH=75
elif [[ "${EXPERIMENT}" == "resub" ]]
then
    CUTOFF=5
    RLENGTH=75
else
    echo "Specify a valid experiment: railGEU, resub"
fi


# Construct shell files
for chrnum in 22 21 Y 20 19 18 17 16 15 14 13 12 11 10 9 8 X 7 6 5 4 3 2 1
do
    chr="chr${chrnum}"
    sname="${SHORT}-${chr}"
    echo "Creating script ${sname}"

    cat > ${ROOTDIR}/.${sname}.sh <<EOF
#!/bin/bash	
#$ -cwd
#$ -m e
#$ -l mem_free=75G,h_vmem=180G,h_fsize=40G
#$ -N ${sname}
#$ -hold_jid fullCov-${EXPERIMENT}

echo "**** Job starts ****"
date

# Make logs directory
mkdir -p ${WDIR}/logs

# Load coverage & get region matrix
cd ${WDIR}
module load R/3.1.x
R -e "maindir <- '${MAINDIR}'; cutoff <- ${CUTOFF}; readLen <- ${RLENGTH}; chr <- '${chr}'; source('${ROOTDIR}/step6-regionMatrix.R')"

## Move log files into the logs directory
mv ${ROOTDIR}/${sname}.* ${WDIR}/logs/

echo "**** Job ends ****"
date
EOF

    call="qsub .${sname}.sh"
    echo $call
    $call
done
