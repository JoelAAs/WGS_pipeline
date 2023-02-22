# Pipeline for Association analysis of genotyped data #

This pipeline was developed as part of the M_Wadelius project.

### Pipeline structure ###

This workflow consists of all the association tests that is to be run on the data. This part is supposed to be run after pipeline_part1 and pipeline_part2, on QCd data.

### Software versions used ###

* bioinfo-tools
* snakemake/4.5.0
* python/3.6.0
* R/3.5.0
* plink/1.90b4.9

### Grouping of case/controls ###

**ADR categories**  
Cases are all samples that has a specific ADR type, and controls are all the rest of the samples. This means controls include all samples with other ADRs, plus population controls were we don't know if they have this particular ADR or not. However, most ADRs are very rare, allowing the assumption that almost none of the controls have the ADR.

**Drug categories**  


**Drug- and ADR categories**  
Cases are all sample that has a specific ADR type on a specific drug type, and controls are all samples other sample on that specific drug (even though they have a different ADR type), plus controls.



### Tests performed by the pipeline ###

The following tests are included in the pipeline:  

__Single variant association analysis on GATK output (SNPs and small Indels)__  
Association analysis in plink using logistic regression with correction for 10 PCs. One test run for each of the inputs in _assoc\_comparisons_, where output is binary (case/control).  

Output from the pipeline:
- run_folder/plink/SNP/plots/{assoc}.png (manhattan plot, one for each test)  
- run_folder/plink/SNP/{assoc}.logistic.assoc (assoc results, one for each test)  

__CYP2D6 metabolizers__  
Association analysis using glm in R, corrected for 10 PCs. Metabolizer status determined using Stargazer, and coded: PM: 0, IM: 1, NM: 2, UM: 3. The model used was:  
_glm(y~metaboliser+PC1+PC2+....+PC10, family = 'binomial')_  
where y was the test specified in the assoc file.  

Output from the pipeline:
- run_folder/CYP2D6/summary_CYP2D6_glm.txt (glm tests on all test performed)  
- run_folder/CYP2D6/R_plots/{assoc}.png (AUC plots for all tests performed)  

### Instructions for running the pipelines ###

**Snakefile_assoc**

This workflow builds on the output from pipeline_part1 and pipeline_part2, and requires the following files to work (if previous workflows have been run this should all be automatic):

* /proj/sens2018106/nobackup/analysis/pipeline_part1/samplefile.txt  
* /proj/sens2018106/nobackup/analysis/pipeline_part1/run_folder/QC/king/merged_kingpc.ped  
* /proj/sens2018106/nobackup/analysis/pipeline_part2/output_data/merged_annovar.hg19_multianno_filtered.vcf.gz  
* /proj/sens2018106/nobackup/analysis/pipeline_part1/run_folder/stargazer/results/cyp2d6_haplotypes.csv  

In addition to these files there is one more file needed:  

* /proj/sens2018106/nobackup/analysis/pipeline_part3/assoc_comparisons  

This file specifies which tests to run association analysis for. To decide on which comparisons to analyse, run this before running the pipeline:  

__python scripts/create_phenofile.py --samplefile /proj/sens2018106/nobackup/analysis/pipeline_part1/samplefile.txt --out phenofile__  

This will produce a file called __phenofile_counts.table__, that contains the numbers for each case/control comparison from the sample file. Decide on which tests to run, and copy over those to the file _assoc\_comparisons_.  


After preparing the _assoc\_comparisons_ file, run the pipeline with:

__sbatch run_pipeline_assoc.sh__
