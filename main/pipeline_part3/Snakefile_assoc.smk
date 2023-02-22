configfile: "config_assoc.yml"


import sys
#sys.path.insert(0, "/proj/sens2018106/softwares/python_packages")
import numpy as np
import re
import os
import pandas as pd

#######################################################################
#                           FUNCTIONS                                 #
#######################################################################

# removes stupid R formatting
def fix_outnames(name):
  new = re.sub('[~/-]', '.', name)

  return new


wildcard_constraints:
    n = "[0-9]+"
#######################################################################
#                           CONFIG                                    #
#######################################################################

# adds which tests to run to config, so no need to put them in config file
fh                   = open("assoc_comparisons", "r")
assoc_list           = [line.strip() for line in fh]
config['assoc_list'] = assoc_list
for file in config['assoc_list']:
  config['output_files'].append("run_folder/plink/SNP/plots/"+file+".png")
#  config['output_files'].append("run_folder/CYP2D6/R_plots/"+fix_outnames(file)+".png")

config["skato_batches"] = 1000
# Get burden test filters sets ## NOTE I give up No clean way to have multiple without d
config["burden_filters"] = os.listdir("input_files/burden_filters")
for burden_set in config["burden_filters"]:
    burden_set = burden_set.replace(".yml", "")
    config["output_files"].append("run_folder/burden/Results/burden_{filter}_results.tsv".format(filter=burden_set))

# Get prevalence
#config["prevalence"] = {}
#with open(config["prevalence_file"], "r") as f:
#    for line in f:
#        association, prevalence_val = line.split("\t")
#        config["prevalence"].update({association: float(prevalence_val)})


######################################################################
#                            INCLUDE                                 #
######################################################################


include: "src/skato.smk"

######################################################################
#                            RULES                                   #
######################################################################

rule all:
    input:
        config["output_files"],
        expand("run_folder/plink/SNP/plots/{assoc}.png", assoc = config['assoc_list']),
        expand("run_folder/freq/{assoc}_{case}.afreq", assoc = config['assoc_list'], case = [1,2])


rule make_pheno_file:
    input:
        "/proj/sens2018106/nobackup/analysis/pipeline_part1/samplefile.txt"
    output:
        "phenofile"
    shell:
        """
        python scripts/create_phenofile.py \
        --samplefile {input} \
        --out {output}
        """

rule make_covar_file:
    params:
        failing_hwe = config["hwe_fail"]
    input:
        vcf = config["filtered_vcf"],
        phenofile = "phenofile"
    output:
        "run_folder/PCA/{assoc}.eigenvec"
    shell:
        """
        plink2 \
          --vcf {input.vcf} \
          --pca \
          --remove input_files/exclude_samples.txt \
          --exclude {params.failing_hwe} \
          --pheno {input.phenofile} \
          --pheno-name {wildcards.assoc} \
          --out run_folder/PCA/{wildcards.assoc} \
          --require-pheno \
          --threads 8
        """

rule get_maf:
    input:
        vcf = config["filtered_vcf"],
        phenofile = "phenofile"
    output:
        freq = "run_folder/freq/{assoc}_{case}.afreq" 
    shell:
        """
        plink2 \
            --vcf {input.vcf} \
            --pheno {input.phenofile} \
            --pheno-name {wildcards.assoc} \
            --freq \
            --keep-if {wildcards.assoc} == {wildcards.case} \
            --set-missing-var-ids @_#:\$r:\$a \
            --new-id-max-allele-len 500 \
            --out run_folder/freq/{wildcards.assoc}_{wildcards.case}
        """


rule get_maf_tot:
    input:
        vcf = config["filtered_vcf"],
        phenofile = "phenofile"
    output:
        freq = "run_folder/freq/{assoc}.afreq"
    shell:
        """
        plink2 \
            --vcf {input.vcf} \
            --pheno {input.phenofile} \
            --pheno-name {wildcards.assoc} \
            --freq \
            --set-missing-var-ids @_#:\$r:\$a \
            --new-id-max-allele-len 500 \
            --out run_folder/freq/{wildcards.assoc}
        """


rule run_plink_SNP_assoc:
    params:
        maf=0.1,
        failing_hwe=config["hwe_fail"]
    input:
        vcf = config["filtered_vcf"],
        covar = "run_folder/PCA/{assoc}.eigenvec",
        phenofile = "phenofile"
    output:
        "run_folder/plink/SNP/{assoc}.{assoc}.glm.logistic.hybrid"
    shell:
        """
        plink2 \
          --vcf {input.vcf} \
          --glm firth-fallback \
          --covar {input.covar} \
          --covar-name PC1 PC2 PC3 PC4 \
          --set-missing-var-ids @_#:\$r:\$a \
          --new-id-max-allele-len 500 \
          --pheno {input.phenofile} \
          --pheno-name {wildcards.assoc} \
          --remove input_files/exclude_samples.txt \
          --exclude {params.failing_hwe} \
          --covar-variance-standardize \
          --threads 10 \
          --out run_folder/plink/SNP/{wildcards.assoc}
    """

rule manhattan_polt:
    input:
        assoc_res = "run_folder/plink/SNP/{assoc}.{assoc}.glm.logistic.hybrid"
    output:
        "run_folder/plink/SNP/plots/{assoc}.png"
    shell:
        """
         python3 scripts/manhattan_plot.py \
          --assoc {input} \
          --out {output}
        """

rule run_glm_CYP2D6:
    input:
        cypfile    = "/proj/sens2018106/nobackup/analysis/pipeline_part1/run_folder/stargazer/results/cyp2d6_haplotypes.csv",
        pcfile     = "/proj/sens2018106/nobackup/analysis/pipeline_part1/run_folder/QC/king/merged_kingpc.ped",
        samplefile = "/proj/sens2018106/nobackup/analysis/pipeline_part1/samplefile.txt",
        assocfile  = "/proj/sens2018106/nobackup/analysis/pipeline_part3/assoc_comparisons",
        phenofile  = "/proj/sens2018106/nobackup/analysis/pipeline_part3/phenofile"
    output:
        "run_folder/CYP2D6/summary_CYP2D6_glm.txt",
        expand("run_folder/CYP2D6/R_plots/{assoc}.png", assoc = [fix_outnames(name) for name in config['assoc_list']])
    shell:
        """
        Rscript --vanilla scripts/analysis_logistic_MW.R {input.cypfile} {input.pcfile} {input.samplefile} {input.assocfile} {input.phenofile} {output[0]} "run_folder/CYP2D6/R_plots/"
        """


############################ GCTA Rules ############################


#rule REML_group:
#    params:
#        gcta = config["gcta"],
#        maf = 0.01
#    input:
#        pheno = "run_folder/GCTA/tests/{assoc}/{assoc}.pheno",
#        multi_grm = "run_folder/GCTA/tests/{assoc}/{assoc}_multi_grm.txt",
#        prevalence = lambda wc: config["prevalence"][wc.assoc]
#    output:
#        "run_folder/GCTA/results/{assoc}/{assoc}.hsq"
#    shell:
#        """
#        {parms.gcta} --reml \
#            --mgrm {input.multi_grm} \
#            --pheno {input.pheno} \
#            --prevalence {input.prevalence} \
#            --maf {params.maf}
#            --out {wildcards.assoc}
#        
#        mv $(echo {output} | sed s"/results/tests/"g) {output}
#        """
#
#rule group_grm_on_ld:
#    params:
#        group_grm_script = config["group_grm"]
#    input:
#        "run_folder/GCTA/tests/{assoc}/{assoc}.score.ld"
#    output:
#        "run_folder/GCTA/tests/{assoc}/{assoc}_multi_grm.txt"
#    shell:
#        """
#        rscript {params.group_grm_script} {input} {wildcards.assoc}
#        """
#
#rule calculate_ld_grm:
#    params:
#        gcta = config["gcta"],
#        segment_kb = 200
#    input:
#        "run_folder/GCTA/tests/{assoc}/{assoc}.bed"
#    output:
#        "run_folder/GCTA/tests/{assoc}/{assoc}.score.ld"
#    shell:
#         """
#         {params.gcta} \
#            --bfile $(echo {input} | sed s"/.bam//"g \
#            --ld-score-region {params.segment_kb} \
#            --out $(echo {input} | sed s"/.bam//"g
#         """
#
#rule generate_bed_and_pheno:
#    input:
#        vcf = "/proj/sens2018106/nobackup/analysis/pipeline_part2/output_data/merged_annovar.hg19_multianno_filtered.vcf.gz",
#        phenofile = "phenofile"
#    output:
#        bed = "run_folder/GCTA/tests/{assoc}/{assoc}.bed",
#        pheno = "run_folder/GCTA/tests/{assoc}/{assoc}.pheno"
#    run:
#        shell('plink --vcf '+ input.vcf +' --no-sex-allowed --out ' + re.sub("\.bed$", "", output.bed) + ' --make-bed')
#        gcta_phenofile = pd.read_table(input.phenofile, sep="\t")
#        gcta_phenofile.to_csv(output.pheno, columns=["FID", "IID", wildcards.assoc], sep="\t", index=False)
