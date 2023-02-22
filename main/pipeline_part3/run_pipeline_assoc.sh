#! /bin/bash -l
#SBATCH -A sens2018106
#SBATCH -p core
#SBATCH -n 1
#SBATCH -J Assoc_pipe
#SBATCH -t 00-05:00:00

module load bioinfo-tools
module load snakemake/4.5.0
module load plink/1.90b4.9
module load python/3.6.0
module load R/3.5.0


snakemake -r -p \
--cluster-config cluster.json \
--cluster "sbatch -A {cluster.account} -p {cluster.partition} -n {cluster.n} -t {cluster.time} -o slurm-logs/slurm-%j.out" \
--jobs 150 \
--snakefile Snakefile_assoc


#snakemake --snakefile Snakefile_rest --dag | dot -Tpng > pipeline_dag_rest.png
#snakemake --snakefile Snakefile_rest --rulegraph | dot -Tpng > pipeline_rules_rest.png

