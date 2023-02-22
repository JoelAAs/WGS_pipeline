#! /bin/bash -l
#SBATCH -A sens2018106
#SBATCH -p core
#SBATCH -n 1
#SBATCH -J Rest_pipe
#SBATCH -t 00-05:00:00

module load bioinfo-tools
module load GATK/3.8-0
module load snakemake/4.5.0
module load htslib/1.8
module load samtools/1.6
module load vcftools/0.1.15
module load bcftools/1.6
module load snpEff/4.3t
module load python3
module load vep/92
module load perl_modules/5.24.1
module load MariaDB/10.1.29
module load annovar/2017.07.16
module load R/3.5.0
module load plink/1.90b4.9


snakemake -r -p \
--cluster-config cluster.json \
--cluster "sbatch -A {cluster.account} -p {cluster.partition} -n {cluster.n} -t {cluster.time} -o slurm-logs/slurm-%j.out" \
--jobs 150 \
--snakefile Snakefile_rest  


#snakemake --snakefile Snakefile_rest --dag | dot -Tpng > pipeline_dag_rest.png
#snakemake --snakefile Snakefile_rest --rulegraph | dot -Tpng > pipeline_rules_rest.png

