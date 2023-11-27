#!/bin/bash
# directories with code

#example call: bash bash/cells/submit_psm.sh 10 20 1.0 0.01 -0.01 0.002 0.002 3000 pi_ohern,day,scavenge 0-12:00:00 1 1
cellsdir=~/dpm
srcdir=$cellsdir/src
maindir=$cellsdir/main/cell

# directory for all output for cell simulations
outputdir=/gpfs/gibbs/pi/ohern/at965/dpm

# directory for simulations specific to psm
simtypedir=$outputdir/psm

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
phi=$4
kl=$5
ka=$6
kb=$7
att=$8
att2=$9
t_maxwell="${10}"
v0="${11}"
t_abp="${12}"
duration="${13}"
partition="${14}"
time="${15}"
numRuns="${16}"
startSeed="${17}"

numSeedsPerRun=1
let numSeeds=$numSeedsPerRun*$numRuns
let endSeed=$startSeed+$numSeeds-1

# name strings
basestr=psm_calA0"$calA0"_phi"$phi"_tm"$t_maxwell"_v0"$v0"_t_abp"$t_abp"_kl"$kl"_ka"$ka"_"kb"$kb
# name of file (not unique, but is unique within its subdirectory)
file_basename=_N"$NCELLS"_dur"$duration"_att"$att"_att2"$att2"_start"$startSeed"_end"$endSeed"
# name of run (must be unique, distinguishes from other runs on slurm)
runstr="$basestr""$file_basename"

# make directory specific for this simulation
simdatadir=$simtypedir/$basestr
mkdir -p $simdatadir

# write input parameters to a configuration file for organization
configFile=$simdatadir/"$file_basename"_config.txt

# compile into binary
binf=bin/"$runstr".o
mainf=$maindir/psm2D.cpp

echo Running psm simulations with parameters: > $configFile
echo NCELLS = "$NCELLS" >> $configFile
echo NV = "$NV" >> $configFile
echo calA0 = "$calA0" >> $configFile
echo phi = "$phi" >> $configFile
echo kl = "$kl" >> $configFile
echo ka = "$ka" >> $configFile
echo kb = "$kb" >> $configFile
echo att = "$att" >> $configFile
echo att2 = "$att2" >> $configFile
echo t_maxwell = "$t_maxwell" >> $configFile
echo v0 = "$v0" >> $configFile
echo t_abp = "$t_abp" >> $configFile
echo k_ecm = "$k_ecm" >> $configFile
echo k_off = "$k_off" >> $configFile
# run compiler
rm -f $binf
g++ --std=c++11 -O3 -I "$srcdir" "$mainf" "$srcdir"/dpm.cpp "$srcdir"/cell.cpp -o $binf
echo compiling with : g++ --std=c++11 -O3 -I "$srcdir" "$mainf" "$srcdir"/dpm.cpp "$srcdir"/cell.cpp -o $binf

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
        filestr="$file_basename"_sd"$seed"

        # create output files
        outFileStem=$simdatadir/$filestr

        # append to runString
        runString="$runString ; ./$binf $NCELLS $NV $calA0 $phi $kl $ka $kb $att $att2 $t_maxwell $v0 $t_abp $seed $duration $outFileStem"
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
echo \#SBATCH --requeue
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
