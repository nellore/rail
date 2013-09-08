#./flux-simulator/bin/flux-simulator -p parameters/sailfish.par
#SIM_OUT=simulation_output
python $TORNADO_HOME/src/simulate/flux.py \
    --output-dir=$SIM_HOME \
    --flux-path=$FLUX_HOME \
    --gtf-file=$PWD/genes.fixed.gtf \
    --chromosomes=$PWD/Chromosomes \
    --num-reads=1500 \
    --num-molecules=5000 \
    --pcr-probability=1 \
    --num-samples=2
    
#appends the file name in front of every read name
#ex: sh format_sim.sh small_test
echo "Reformatting simuated data"
for FASTQ in $SIM_HOME/*.fastq
do
    name=`basename $FASTQ .fastq`
    echo $name
    sed -i '1~4 s/\@/\@'$name'./' $FASTQ
done
