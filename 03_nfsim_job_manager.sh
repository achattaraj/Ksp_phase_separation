#!/bin/bash
#SBATCH --job-name=FTC_60uM
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --array=1-100
#SBATCH --partition=general
#SBATCH --qos=specialuse
#SBATCH --mail-type=END
#SBATCH --mem=1G
#SBATCH -o /home/CAM/achattaraj/NFSim/jobfiles/%x_run%a.out


inpath="/home/CAM/achattaraj/NFSim"

sysname="4v_4v_FTC_60uM"

outpath="$inpath/$sysname"
stdpath="$inpath/$sysname/stdout"

mkdir -p $outpath
mkdir -p $stdpath

t_end=0.5
numSteps=100
utl=2

module load BioNetGen/2.5.0

NFsim "-xml" "$inpath/$sysname.xml" "-sim" "$t_end" "-v" "-seed" "$SLURM_ARRAY_TASK_ID" "-o" "$outpath/Run_$SLURM_ARRAY_TASK_ID.gdat" "-ss" "$outpath/Run_$SLURM_ARRAY_TASK_ID.species" "-utl" "$utl" "-cb" "-oSteps" "$numSteps" 1>"$stdpath/Run_$SLURM_ARRAY_TASK_ID.out" 2>"$stdpath/Run_$SLURM_ARRAY_TASK_ID.err"
