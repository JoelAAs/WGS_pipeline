#!/usr/bin/env python3
import sys

sys.path.insert(0, "/proj/sens2018106/softwares/python_packages")
import tabix
import os
import pandas as pd
import re
import yaml
import multiprocessing as mp


def add_genotype(sample):
    """ Sum genotypes """
    if sample[0] is not ".":
        return int(sample[0]) + int(sample[2])
    else:
        return None


def get_genotype_sum(genotypes):
    """ Get genotypes and sum allels for matrix and maf"""
    allel_list = list(map(add_genotype, genotypes))
    allele_count = sum(filter(None, allel_list))
    return allel_list, allele_count


def filter_annotation_data(ann_str, fields):
    def _get_value(annotation, field_name):
        pattern = "(?<=" + field_name + "=)[\w.]+(?=;)"
        re_match = re.search(pattern, annotation)
        if re_match:
            return re_match.group()
        return "."

    values = [_get_value(ann_str, curr_field) for curr_field in fields]
    return values


def write_genotype_lines(filename_bgz, curr_chrom, from_tp, to_tp, annotation_fields, gmat_file, var_id_file):
    tabix_vcf = tabix.open(filename_bgz)
    variant_lines = tabix_vcf.query(curr_chrom, from_tp, to_tp)
    var_pk = 0
    for variant_line in variant_lines:
        annotations = filter_annotation_data(variant_line[7], annotation_fields)
        try:
            matrix_line, allele_count = get_genotype_sum(variant_line[9:])
        except ValueError:
            raise ValueError(variant_line[9:])

        # Format
        matrix_line = "\t".join(map(lambda x: str(x) if x is not None else "", matrix_line))
        var_id_file.write(
            "\t".join([str(var_pk)] +
                      variant_line[:5] +  # CHROM POS ID REF ALT
                      [str(allele_count)] +
                      annotations) + "\n")

        gmat_file.write(str(var_pk) +
                        "\t" + matrix_line + "\n")
        var_pk = var_pk + 1


def write_gmat_var_pos(filename_vcf, burden_wd, filter_fields, samples, chrom, from_bp, to_bp, gene_name, filter_name):
    # Collection job for multiprocessing
    gmat_file = open(f"{burden_wd}{filter_name}/variant_matrix/gmats/{gene_name}_gmat.tsv", "w+")
    gmat_file.write("PK\t" + samples + "\n")
    var_id_file = open(f"{burden_wd}{filter_name}/variant_matrix/vars/{gene_name}_var.tsv", "w+")
    var_id_file.write("PK\tCHROM\tPOS\tID\tREF\tALT\tCOUNT\t" + "\t".join(filter_fields) + "\n")
    write_genotype_lines(filename_vcf, chrom, from_bp, to_bp, filter_fields, gmat_file, var_id_file)
    gmat_file.close()
    var_id_file.close()


if __name__ == '__main__':
    # Arguments
    args = sys.argv[1:]
    filename_bgz = args[0]
    gene_set = pd.read_csv(args[1], sep="\t", dtype={"CHROM": str})
    gene_set = gene_set[~gene_set.gene_name.duplicated()]
    filter_config_filename = args[2]
    with open(filter_config_filename, "r") as ymlfile:
        config = yaml.load(ymlfile)

    filter_name = os.path.basename(filter_config_filename).replace(".yml", "")

    assert config["burden_filters"] is not None, "No 'burden_filters' in config"

    # What fields are of interest
    filter_fields = []

    for filter_type in config["burden_filters"]:
        part_filters = config["burden_filters"][filter_type]
        if filter_type == "all":
            pass
        else:
            if part_filters:
                filter_fields = filter_fields + list(part_filters.keys())
    filter_fields = filter_fields + config["burden_annotation"]
    filter_fields = list(set(filter_fields))

    # Check if TBI exists in same folder
    if not os.path.isfile(filename_bgz + ".tbi"):
        print(filename_bgz + ".tbi: Not found.\n")
        os.system("tabix -p vcf " + filename_bgz)

    cwd = os.getcwd()

    os.system("mkdir -p " + cwd + f"/run_folder/burden/{filter_name}/variant_matrix/gmats")
    os.system("mkdir -p " + cwd + f"/run_folder/burden/{filter_name}/variant_matrix/vars")
    burden_wd = "run_folder/burden/"

    # Get samples
    samples = os.popen("vcf-query -l " + filename_bgz).read()
    samples = samples.strip().replace("\n", "\t")
    print("Found {} samples in VCF".format(len(samples.split("\t"))))

    allowed_chroms = [str(i) for i in range(1, 24)]
    allowed_chroms = allowed_chroms + ["X", "Y"]
    gene_set = gene_set[gene_set["CHROM"].isin(allowed_chroms)]
    gene_zip = zip(gene_set.CHROM, gene_set.START, gene_set.END, gene_set.gene_name)

    jobs = []
    my_pool = mp.Pool(4)

    for chrom, from_bp, to_bp, gene_name in gene_zip:
        job = my_pool.apply_async(
            write_gmat_var_pos,
            args=(filename_bgz, burden_wd, filter_fields, samples, chrom, from_bp, to_bp, gene_name, filter_name))
        jobs.append(job)

    for job in jobs:
        print(job.get())

    my_pool.close()
