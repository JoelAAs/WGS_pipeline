# Script for filtering data #

This script is used to filter the files produced by the other pipelines, before being used for association analysis.  

Run it as it is or step by step.  

By default samples to exclude and snps to exclude are taken from the QC results produced by king. What is not in there is PCA outliers, those have to be added manually to the samples_to_exclude.txt file. This is easily done when running the script step by step.  

This script uses the output from the Annovar annotated VCF file only.  

Before running, softlink all resulting VCF files from other pipelines to the folder input_data.  

The VCF file is unzipped and rezipped when filtered, just because the type of zipping is not really compatible with the filtering.  


__Run with:__  
./filter_vcf.sh 
