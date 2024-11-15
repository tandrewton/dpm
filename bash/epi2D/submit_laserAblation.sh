#!/bin/bash
# directories with code

#example call: bash bash/epi2D/submit_laserAblation.sh 20 20 4 1.10 0.92 0.925 1.0 1.0 0.01 0.1 0.01 0.0 1.0 2.0 1.0 3.0 10.0 0 0 1 200 pi_ohern,day,scavenge 0-12:00:00 1 1
#weird bug with configFile
cellsdir=~/dpm
srcdir=$cellsdir/src
maindir=$cellsdir/main/epi2D

# directory for all output for cell simulations
outputdir=/gpfs/gibbs/pi/ohern/at965/dpm

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
ka=$8
kb=$9
att="${10}"
strainRate_ps="${11}"
deltaSq="${12}"
k_ps="${13}"
k_lp="${14}"
tau_lp="${15}"
d_flag="${16}"
t_stress="${17}"
taus="${18}"
smooth="${19}"
duration="${20}"
partition="${21}"
time="${22}"
numRuns="${23}"
startSeed="${24}"

numSeedsPerRun=1
let numSeeds=$numSeedsPerRun*$numRuns
let endSeed=$startSeed+$numSeeds-1

# name strings 
basestr=ablate_A0"$calA0"_t_stress"$t_stress"k_l"$kl"_k_a"$ka"_k_b"$kb"_w_ps"$strainRate_ps"_dsq"$deltaSq"_k_ps"$k_ps"_k_lp"$k_lp"_d_flag"$d_flag"_taus"$taus"_sm"$smooth"

# name of file (not unique, but is unique within its subdirectory)
file_basename=_N"$NCELLS"_Dur"$duration"_att"$att"_sd"$startSeed"_sd"$endSeed"

runstr="$basestr""$file_basename"

# make directory specific for this simulation
simdatadir=$simtypedir/$basestr
mkdir -p $simdatadir

# write input parameters to a configuration file for organization
configFile=$simdatadir/"$file_basename"_config.txt

# compile into binary
binf=bin/"$runstr".o
mainf=$maindir/laserAblation.cpp

echo Running laserAblation simulations with parameters: >> $configFile
echo NCELLS = "$NCELLS" >> $configFile
echo NV = "$NV" >> $configFile
echo ndelete = "$ndelete" >> $configFile
echo calA0 = "$calA0" >> $configFile
echo phiMin = "$phiMin" >> $configFile
echo phiMax = "$phiMax" >> $configFile
echo kl = "$kl" >> $configFile
echo ka = "$ka" >> $configFile
echo kb = "$kb" >> $configFile
echo att = "$att" >> $configFile
echo strainRate_ps = "$strainRate_ps" >> $configFile
echo deltaSq = $deltaSq >> $configFile
echo k_ps = $k_ps >> $configFile
echo k_lp = $k_lp >> $configFile
echo tau_lp = $tau_lp >> $configFile
echo d_flag = $d_flag >> $configFile
echo t_stress = $t_stress >> $configFile
echo taus = "$taus" >> $configFile
echo smooth = "$smooth" >> $configFile
echo duration = "$duration" >> $configFile


# run compiler
rm -f $binf
g++ --std=c++11 -O3 -I "$srcdir" "$mainf" "$srcdir"/dpm.cpp "$srcdir"/epi2D.cpp -o $binf
echo compiling with : g++ --std=c++11 -O3 -I "$srcdir" "$mainf"  "$srcdir"/dpm.cpp "$srcdir"/epi2D.cpp -o $binf

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
        runString="$runString ; ./$binf $NCELLS $NV $ndelete $calA0 $phiMin $phiMax $kl $ka $kb $att $strainRate_ps $deltaSq $k_ps $k_lp $tau_lp $d_flag $t_stress $taus $smooth $seed $duration $outFileStem"
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
# 3. number of cells to delete
# 3. calA0
# 4. phiMin
# 5. phiMax
# 6. kl
# 7. ka
# 8. att
# 9. segment shrinking rate for pursestring
# 10. max yield length for purse string springs
# 11. spring constant purse string virtual particle to wound vertex
# 12. spring constant for flag to nearest vertex on wound edge for crawling
# 13. protrusion time constant (controls stochastic lifetime of a protrusion)
# 14. protrusion distance from cell edge
# 16. rotational diffusion constant for cell's protrusion director
# 17. whether cells have CIL or not
# 18. number of timesteps
# 19. partition
# 20. time
# 21. number of runs (number of array entries, i.e. arraynum)
# 22. start seed (end seed determined by number of runs)
