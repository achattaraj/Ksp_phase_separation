#!/bin/bash
#SBATCH --job-name=RefSys
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --array=0-49
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=1G
#SBATCH -o /home/CAM/achattaraj/Springsalad/stdout/%x_%A_Run%a.out
#SBATCH -e /home/CAM/achattaraj/Springsalad/stdout/%x_%A_Run%a.err

module load java/1.8.0_77

java -Xms64m -Xmx1024m -jar /home/CAM/achattaraj/Springsalad/jar/LangevinNoVis01.jar "/home/CAM/achattaraj/Springsalad/sims/A4_B4_reference_system_SIM.txt" "$SLURM_ARRAY_TASK_ID" 2> /home/CAM/achattaraj/Springsalad/sims/out

