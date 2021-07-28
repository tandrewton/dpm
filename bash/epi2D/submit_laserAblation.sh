#!/bin/bash
# directories with code

#example call: bash bash/epi2D/submit_laserAblation.sh 24 24 5 1.08 0.7 0.9 1.0 0.5 0.4 0.5 0.0 1 1000 pi_ohern,day,scavenge 0-12:00:00 1 1
cellsdir=~/dpm
srcdir=$cellsdir/src
maindir=$cellsdir/main/epi2D

# directory for all output for cell simulations
outputdir=/gpfs/loomis/project/fas/ohern/at965/dpm

# directory for simulations specific to laserAblation
simtypedir=$outputdir/ablate

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
ndelete=$3
calA0=$4
phiMin=$5
phiMax=$6
kl=$7
att=$8
v0=$9
B="${10}"
Dr0="${11}"
boolCIL="${12}"
NT="${13}"
partition="${14}"
time="${15}"
numRuns="${16}"
startSeed="${17}"

numSeedsPerRun=1
let numSeeds=$numSeedsPerRun*$numRuns
let endSeed=$startSeed+$numSeeds-1

# name strings
basestr=ablate_N"$NCELLS"_NV"$NV"_ndel"$ndelete"_calA0"$calA0"_kl"$kl"_att"$att"_v0"$v0"_B"$B"_Dr0"$Dr0"_CIL"$boolCIL"_NT"$NT"
runstr="$basestr"_startseed"$startSeed"_endseed"$endSeed"

# make directory specific for this simulation
simdatadir=$simtypedir/$basestr
mkdir -p $simdatadir

# compile into binary using packing.h
binf=bin/"$runstr".o
mainf=$maindir/laserAblation.cpp

echo Running laserAblation simulations with parameters:
echo NCELLS = "$NCELLS"
echo NV = "$NV"
echo ndelete = "$ndelete"
echo calA0 = "$calA0"
echo phiMin = "$phiMin"
echo phiMax = "$phiMax"
echo kl = "$kl"
echo att = "$att"
echo v0 = "$v0"
echo B = "$B"
echo Dr0 = "$Dr0"
echo boolCIL = "$boolCIl"
echo NT = "$NT"


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
        filestr="$basestr"_seed"$seed"

        # create output files
        posf=$simdatadir/$filestr.pos
        energyf=$simdatadir/$filestr.energy
        stressf=$simdatadir/$filestr.stress

        # append to runString
        runString="$runString ; ./$binf $NCELLS $NV $ndelete $calA0 $phiMin $phiMax $kl $att $v0 $B $Dr0 $boolCIL $runseed $NT $posf $energyf $stressf"
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
# 3. number of cells to delete
# 3. calA0
# 4. phiMin
# 5. Ptol
# 6. kl
# 7. kb
# 8. att
# 9. active velocity scale
# 10. rotational diffusion constant
# 10. whether cells have CIL or not
# 10. number of timesteps
# 11. partition
# 12. time
# 13. number of runs (number of array entries, i.e. arraynum)
# 14. start seed (end seed determined by number of runs)
