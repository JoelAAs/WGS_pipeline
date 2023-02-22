rule all:
    input:
        "output_data/merged_annovar.hg19_multianno_filtered.vcf.gz"
#        "all_swegene_pg"

rule getidx:
    input:
        "input_data/merged_annovar.hg19_multianno.vcf.gz"
    output:
        "input_data/merged_annovar.hg19_multianno.vcf.gz.tbi"
    shell:
         """
         tabix {input}
         """

rule set_missing_vars:
    input:
        "input_data/merged_annovar.hg19_multianno.vcf.gz",
        "input_data/merged_annovar.hg19_multianno.vcf.gz.tbi"
    output:
        "run_folder/bed/hg19_multianno_varset.bed"
    shell:
         """
         plink2 \
            --vcf {input[0]} \
            --set-missing-var-ids @_#_\$r:\$a \
            --new-id-max-allele-len 300 missing \
            --out run_folder/bed/hg19_multianno_varset \
            --make-bed \
            --threads 10
         """


rule get_faulty_variants:
    params:
         swegen_samples="qc_files/swedegene.txt"
    input:
         bed = "run_folder/bed/hg19_multianno_varset.bed"
    output:
          "run_folder/failed_variants/snps_to_remove_hardy_cr95"
    run:
        plink_input=input.bed.replace(".bed", "")
        shell(f"plink2 --bfile {plink_input}" +
              " --missing --out run_folder/failed_variants/flag_missing")
        shell("cat run_folder/failed_variants/flag_missing.vmiss |" +
              " awk '{{if($5>0.05)print$0}}' | sort -k5 >" +
              " run_folder/failed_variants/snps_call_rate_below_0.95")
        shell(f"plink2 --bfile {plink_input} --keep {params.swegen_samples} --hardy" +
              " --out run_folder/failed_variants/hwe_failed")
        shell("cat run_folder/failed_variants/hwe_failed.hardy |" +
              " awk '{{if($10<1e-8)print$0}}' |" +
              " sort -k10g >" +
              " run_folder/failed_variants/snps_failing_hardy")

        shell("cut -f2 run_folder/failed_variants/snps_call_rate_below_0.95 > " +
              f"{output}")
        shell("cut -f2 run_folder/failed_variants/snps_failing_hardy >> " +
              f"{output}")

        shell(f"sort {output} | uniq > temp")
        shell(f"mv temp {output}")


rule get_all_swegene_pg:
    output:
        "all_swegene_pg"
    shell:
        """
        for file in /dataset/swegen/hg19/BAM/*.bam; do
            samtools view -H "$file" | grep "ID:bwa[[:space:]]" >> all_swegene_pg
            echo "$(basename "$file")"
        done
        """

rule exclude_plink_var:
    input:
        bed = "run_folder/bed/hg19_multianno_varset.bed",
        exclude_snps = "run_folder/failed_variants/snps_to_remove_hardy_cr95"
    output:
        "output_data/hg19_multianno_varset_hwe_call.vcf"
    run:
        output_vcf = output[0].replace(".vcf", "")
	input_bed = input.bed.replace(".bed", "")
        cmd = f"plink2 --bfile {input_bed} " + \
            f"--exclude {input.exclude_snps} " + \
            f"--out {output_vcf} --recode vcf"
        shell(cmd)


rule reformat_plink_output:
    input:
        exclude_snps = "run_folder/failed_variants/snps_to_remove_hardy_cr95"
    output:
        "run_folder/failed_variants/snp_to_remove"
    run:
        with open(output[0], "w") as w:
            with open(input.exclude_snps, "r") as f:
                for line in f:
                    print(line)
                    chrom, pos, ref_alt = line.strip().split("_")
                    ref, alt = ref_alt.split(":")
                    w.write(f"{chrom}:{pos}:{alt}:{ref}\n")

rule unzip_sort_annovar:
    input:
        "input_data/merged_annovar.hg19_multianno.vcf.gz"
    output:        
        temp("output_data/merged_annovar.hg19_multianno_sorted.vcf")
    shell:
        """
        bcftools sort -m 40G -O v --temp-dir output_data/tmp  {input} > {output}
        """


rule main_filter:
    input:
        vcf = "output_data/merged_annovar.hg19_multianno_sorted.vcf",
        snp_list = "run_folder/failed_variants/snp_to_remove",
        sample_list = "qc_files/samples_to_remove.txt"
    output:
        vcf = temp("output_data/merged_annovar.hg19_multianno_filtered.vcf")
    shell:
        """
        python3 scripts/filter_vcf.sys \
            {input.vcf} \
            {input.snp_list} \
            {input.sample_list} \
            {output.vcf}      
        """

rule bgzip:
    input:
        vcf = "output_data/merged_annovar.hg19_multianno_filtered.vcf"
    output:
        vcf = temp("output_data/merged_annovar.hg19_multianno_filtered.vcf.gz")
    shell:
        """
        bgzip {input.vcf}
        tabix {output.vcf}
        """

