#! /bin/bash -l
#SBATCH -A sens2018106
#SBATCH -p core
#SBATCH -n 1
#SBATCH -J Name_run
#SBATCH -t 20:00:00

module load bioinfo-tools
module load GATK/3.8-0
module load snakemake/4.5.0
module load python3


snakemake -r -p \
--cluster-config cluster.json \
--cluster "sbatch -A {cluster.account} -p {cluster.partition} -n {cluster.n} -t {cluster.time} -o slurm-logs/slurm-%j.out" \
--jobs 50 \
--snakefile Snakefile_batching


#snakemake --snakefile Snakefile_batching --dag | dot -Tpng > pipeline_dag.png
#snakemake --snakefile Snakefile_batching --rulegraph | dot -Tpng > pipeline_rules.png

