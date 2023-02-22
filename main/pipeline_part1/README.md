# Pipeline for GATK variant calling #

This pipeline was developed as part of the M_Wadelius project.

### Pipeline structure ###

This workflow consists of a few separate pipelines intended to be run in a specified order, with separate Snakefiles:

* Snakefile_batching
* Snakefile_rest

**Snakefile_batching**  
This pipeline splits samples into batches, and also into genomic intervals, and runs CombineGVCF on all batches. This pipeline can be re-run as soon as new samples have been added, without risking removing any samples already processed.

**Snakefile_rest**  
This pipeline is meant to be run when all samples have been processed through the first pipeline. This pipeline runs GenotypeGVCFs, combines the intervals, runs VQSR for SNPs and Indels, removes failed variants, and splits multiallelic sites. In addition it runs the following QC tests: BamQC, SamtoolsStats, MultiQC, king for various variant metrics, and BcfTools. Finally it does annotation using SnpEff, VEP, and Annovar.

### Software versions used ###

* bioinfo-tools
* GATK/3.8-0
* snakemake/4.5.0
* htslib/1.8
* samtools/1.6
* vcftools/0.1.15
* bcftools/1.6
* snpEff/4.3t
* python3
* vep/92
* perl_modules/5.24.1
* MariaDB/10.1.29
* annovar/2017.07.16
* R/3.5.0
* plink/1.90b4.9


### Instructions for running the pipelines ###

**Snakefile_batching**

Files to be edited before run:
* samplelist.txt: a list of all samples to be processed, in the same format as the example provided here.
* config_batching.yml: a yaml file specifying the input parameters used, plus the names of the g.vcf files to be processed


When running the pipeline for the first time, set batch to 1, and select a good interval size, this cannot be changed later. Run the pipeline with sbatch run_pipeline_batches.sh, which will submit the job to slurm.

When receiving additional samples, add these to the samplefile.txt and config_batching.yml file (no need to remove already processed samples). Change the batch number in the config_batching.yml file and submit the pipeline to slurm. Before running it will first check that the interval size is the same as previous runs, and that the batch number is unique (this prevents it from overwriting any already processed files). The bin size however, can be different between batches without it affecting the results.

**Snakefile_rest**

After running all the samples through the first pipeline, run this pipeline once.

No files have to be edited before running this pipeline (although the samplefile.txt has to still contain all the samples), it will look up all the batch, bin, and interval files created by the first pipeline and run GenotypeGVCFs on all of them. Simply do sbatch run_pipeline_rest to submit it to slurm.

### Output files from the pipelines ###

The following output files are the files of most interest:

* Unannotated VCF file:  
   \- run_folder/gatk/gvcfs/joint/geno/VQSR/filter/merged_VQSR_variants_pass_split.vcf.gz


* SnpEff annotated VCF file:  
  \- run_folder/gatk/gvcfs/joint/geno/VQSR/filter/annotated/merged_snpEff.vcf.gz


* VEP annotated txt file:  
  \- run_folder/gatk/gvcfs/joint/geno/VQSR/filter/annotated/merged_VEP.txt.gz


* Annovar annotated VCF file:  
  \- run_folder/gatk/gvcfs/joint/geno/VQSR/filter/annotated/merged_annovar.hg19_multianno.vcf.gz


* MultiQC report:  
  \- run_folder/QC/multiQC/multiqc_report.html


* QC metrics from king:  
  \- run_folder/QC/king/merged_king_


* BcfTools metrics:  
  \- run_folder/QC/bcftools/merged_bcftools_stats.txt