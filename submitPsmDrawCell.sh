#!/bin/bash
#SBATCH --job-name myjob
#SBATCH --cpus-per-task 4
#SBATCH --mem 18G
#SBATCH -t 8:00:00

module load MATLAB/2023a
# assuming you have your_script.m in the current directory
matlab -batch "drawCellSim("40", "1.0", "0.8", "5.0", "0.1", "0.001", "0.001", "0.1", "10000.0", "0")"