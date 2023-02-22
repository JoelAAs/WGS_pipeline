#! /bin/bash -l
#SBATCH -A sens2018106
#SBATCH -p core
#SBATCH -n 2
#SBATCH -J filter
#SBATCH -t 00-24:00:00


# extract which samples to remove
grep -v IID /proj/sens2018106/nobackup/analysis/pipeline_part1/run_folder/QC/king/merged_king_autoQC_sampletoberemoved.txt \
| cut -f1 > samples_to_remove.txt

# extract which snps to exclude. removes all from king analysis (call rate 95%)
grep -v SNP /proj/sens2018106/nobackup/analysis/pipeline_part1/run_folder/QC/king/merged_king_autoQC_snptoberemoved.txt \
| cut -f1 | sed 's/_/:/g' > snps_to_remove.txt

# unzip file because there's something with the zipping the filter doesn't like
echo "$(date)"
gunzip -c input_data/merged_annovar.hg19_multianno.vcf.gz > output_data/temp.vcf

# filter vcf file
echo "$(date)"
./scripts/filter-vcf \
--vcf output_data/temp.vcf \
--list snps_to_remove.txt \
--sample samples_to_remove.txt \
--exclude \
--remove \
--out output_data/merged_annovar.hg19_multianno_filtered.vcf

# gxip output file
echo "$(date)"
gzip output_data/merged_annovar.hg19_multianno_filtered.vcf

#remove temp file
echo "$(date)"
rm output_data/temp.vcf
