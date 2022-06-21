#!/bin/bash
# directories with code

#example call: bash bash/epi2D/submit_neuralTube.sh 10 20 1.0 0.01 0.002 0.002 1 3000 pi_ohern,day,scavenge 0-12:00:00 1 1
dpmdir=~/dpm
srcdir=$dpmdir/src
maindir=$dpmdir/main/cells

# directory for all output for cell simulations
outputdir=/gpfs/loomis/project/fas/ohern/at965/dpm

# directory for simulations specific to laserAblation
simtypedir=$outputdir/cells

# make directories, unless they already exist
mkdir -p $outputdir
mkdir -p $simtypedir
mkdir -p bin
mkdir -p tasks
mkdir -p slurm
mkdir -p out

# inputs
# NT = simulation time (in tau), time = human time
NCELLS=$1
NV=$2
calA0=$3
att=$4
prate=$5
adhrate=$6
duration=$7
partition=$8
time=$9
numRuns="${10}"
startSeed="${11}"

numSeedsPerRun=1
let numSeeds=$numSeedsPerRun*$numRuns
let endSeed=$startSeed+$numSeeds-1

# name strings
basestr=NT_calA0"$calA0"_prate"$prate$"_adhrate"$adhrate$"
runstr="$basestr"_NCELLS"$NCELLS"_Duration"$duration"_att"$att"_startseed"$startSeed"_endseed"$endSeed"

# make directory specific for this simulation
simdatadir=$simtypedir/$basestr
mkdir -p $simdatadir

# write input parameters to a configuration file for organization
configFile=$simdatadir/config.txt

# compile into binary
binf=bin/"$runstr".o
mainf=$maindir/neuralTube.cpp

echo Running neuralTube simulations with parameters: > configFile
echo NCELLS = "$NCELLS" >> $configFile
echo NV = "$NV" >> $configFile
echo calA0 = "$calA0" >> $configFile
echo att = "$att" >> $configFile
echo prate = "$prate" >> $configFile
echo duration = "$duration" >> $configFile
echo prate = "$prate" >> $configFile
echo adhrate = "$adhrate" >> $adhrate

# run compiler
rm -f $binf
g++ --std=c++11 -O3 -I "$srcdir" "$mainf" "$srcdir"/*.cpp -o $binf
echo compiling with : g++ --std=c++11 -O3 -I "$srcdir" "$mainf" "$srcdir"/*.cpp -o $binf

# check compilation
if [[ ! -f $binf ]]
then
    echo -- binary file does not exist, compilation failed.
    exit 1
fi

# create task file
taskf=tasks/"$runstr".task
rm -f $taskf

# loop over files
let fcount=0

# LOOP OVER FILES.
for seed in `seq $startSeed $numSeedsPerRun $endSeed`; do
    # count files
    let fcount=$fcount+1

    # echo to console
    echo On base seed $seed

    # echo string of numSeedPerRun commands to task file
    runString="cd `pwd`"

    # loop over seeds to go into runString
    let ssMax=$numSeedsPerRun-1

    for ss in `seq 0 $ssMax`; do
        # get seed for actual run
        let runseed=$seed+ss

        # get file str
        filestr="$runstr"_seed"$seed"

        # create output files
        outFileStem=$simdatadir/$filestr

        # append to runString
        runString="$runString ; ./$binf $NCELLS $NV $calA0 $att $prate $adhrate $seed $duration $outFileStem"
    done

    # finish off run string
    runString="$runString ;"

    # echo to task file
    echo "$runString" >> $taskf
done

# test if task file was created
if [[ ! -f "$taskf" ]]
then
    echo task file not created, ending before job submission
    exit 1
fi

# get number of jobs to submit to each array
let arraynum=$fcount
echo -- total number of array runs = $arraynum

# setup slurm files
slurmf=slurm/"$runstr".slurm
job_name="$runstr"
runout=out/"$runstr"-%a.out
rm -f $slurmf

# echo about time
echo -- running time = $time for $partition

echo -- PRINTING SLURM FILE...
echo \#\!/bin/bash >> $slurmf
echo \#SBATCH --cpus-per-task=1 >> $slurmf
echo \#SBATCH --array=1-$arraynum >> $slurmf
echo \#SBATCH -n 1 >> $slurmf
echo \#SBATCH -p $partition >> $slurmf
echo \#SBATCH -J $job_name >> $slurmf
echo \#SBATCH --mail-type=END,FAIL >> $slurmf
echo \#SBATCH --mail-user=andrewtondata@gmail.com >> $slurmf
echo \#SBATCH -o $runout >> $slurmf
echo sed -n \"\$\{SLURM_ARRAY_TASK_ID\}p\" "$taskf" \| /bin/bash >> $slurmf
cat $slurmf

# run sbatch file
echo -- running on slurm in partition $partition
sbatch -t $time $slurmf

# ====================
#       INPUTS
# ====================
# 1. NCELLS
# 2. NV
# 3. calA0
# 4. att
# 19. partition
# 20. time
# 21. number of runs (number of array entries, i.e. arraynum)
# 22. start seed (end seed determined by number of runs)
