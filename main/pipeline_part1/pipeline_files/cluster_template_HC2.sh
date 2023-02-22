#! /bin/bash -l
#SBATCH -A sens2018106
#SBATCH -p core
#SBATCH -n 4
#SBATCH -J HC_2
#SBATCH -t 30:00:00
#SBATCH -o slurm-logs/slurm-%j.out

module load bioinfo-tools
module load GATK/3.8-0
module load snakemake/4.5.0
module load python3
