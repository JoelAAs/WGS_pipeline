from snakemake.io import expand

wildcard_constraints:
  burden_filter = "[a-zA-Z0-9-_]+",
  gene = "[a-zA-Z0-9-_.]+"

def get_gene_set(burden_filter):
    filter_file = f"input_files/burden_filters/{burden_filter}.yml"
    gene_set = None
    with open(filter_file, "r") as f:
         for line in f:
            tmp = re.search('gene_tsv.*\"(.*)\"', line)
            if tmp:
                gene_set = tmp.group(1)
                break
    if not gene_set:
        raise ValueError("no geneset in filter")
    return gene_set

def get_genes(wildcards):
    gene_set = get_gene_set(wildcards.burden_filter)
    if gene_set:
        with open(gene_set, "r") as f:
            lines = f.readlines()
            genes = [l.split("\t")[4] for l in lines if not re.search("(MT|GL*)", l.split("\t")[0])][1:]  # remove header
    else:
        raise ValueError("no geneset in filter")
    return genes

def get_input_gene_output(wildcards):
    genes = get_genes(wildcards)
    burden_filter = wildcards.burden_filter

    checkpoint_output = checkpoints.genotype_matrix.get(burden_filter=wildcards.burden_filter).output

    gmats = expand("run_folder/burden/{burden_filter}/variant_matrix/gmats/{gene}_gmat.tsv",
        gene=genes,burden_filter=burden_filter)
    vars = expand("run_folder/burden/{burden_filter}/variant_matrix/vars/{gene}_var.tsv",
        gene = genes, burden_filter=burden_filter)
    shell(f"mkdir -p run_folder/burden/{burden_filter}")
    with open(f"run_folder/burden/{burden_filter}/all.txt", "w")  as w:
        for g, v in zip(gmats,vars):
            w.write(g + "\n")

    return f"run_folder/burden/{burden_filter}/all.txt"

def get_batches_output(wildcards):
    checkpoint_output = checkpoints.get_batches.get(burden_filter=wildcards.burden_filter).output
    n_genes = len(get_genes(wildcards))
    n_batches = int(n_genes / config["skato_batches"]) + 1

    return expand("run_folder/burden/{burden_filter}/batches/{burden_filter}_batch_results_{n}",
        burden_filter=wildcards.burden_filter, n = range(n_batches))

def get_genes_in_batch(wildcards):
    batch_list = "run_folder/burden/{filter}/batches/{filter}_batch_{n}".format(
        filter=wildcards.burden_filter, n=wildcards.n
    )
    with open(batch_list, "r") as f:
        genes = [re.search("\/([a-zA-Z0-9-_.]*)_gmat.tsv$", g.strip()).group(1) for g in f.readlines()]
    return genes


def get_distance_vars_in_batch(wildcards):
    genes = get_genes_in_batch(wildcards)
    distance_vars = [
        f"run_folder/burden/{wildcards.burden_filter}/variant_matrix/annotated_vars/{wildcards.n}/{gene}_var_distance.tsv"
        for gene in genes
    ]
    checkpoint_output = checkpoints.get_variants_in_genes.get(**wildcards).output
    return distance_vars


rule concat_skato_results:
    input:
        lambda wc: get_batches_output(wc)
    output:
        "run_folder/burden/Results/burden_{burden_filter}_results.tsv"
    run:
        with open(output[0], "w") as w:
            for fname in input:
                with open(fname.strip(),"r") as f:
                    for line in f:
                        with open(line.strip(), "r") as res:
                            for l in res:
                                w.write(l)


rule skato_per_gene:
    """ If we add assoc_comparisons to file change coded path here """
    input:
        batch_list="run_folder/burden/{burden_filter}/batches/{burden_filter}_batch_{n}",
        gmats_in_batch = lambda wc: expand("run_folder/burden/{{burden_filter}}/variant_matrix/gmats/{gene}_gmat.tsv", gene = get_genes_in_batch(wc)),
        anno_var = get_distance_vars_in_batch,
        phenotype = "phenofile",
        burden_filter = "input_files/burden_filters/{burden_filter}.yml",
        chosen_tests ="assoc_comparisons",
        excluded_samples = "input_files/exclude_samples.txt",
        covar_file = expand("run_folder/PCA/{assoc}.eigenvec", assoc = config['assoc_list'])
    output:
        "run_folder/burden/{burden_filter}/batches/{burden_filter}_batch_results_{n}"
    run:
        shell(f"mkdir -p run_folder/burden/{wildcards.burden_filter}/results")
        with open(output[0], "w") as w:
            for distance_file in input.anno_var:
                cmd = f"Rscript scripts/analysis_skato.R \
                        {distance_file} \
                        {input.phenotype} \
                        {input.burden_filter} \
                        {input.chosen_tests} \
                        {input.excluded_samples} \
                        {input.covar_file}"
                shell(cmd)
                gene = re.search("[a-zA-Z0-9]+_var_distance.tsv", distance_file).group(0)
                gene = gene.replace("_var_distance.tsv", "_result.tsv")
                w.write(f"run_folder/burden/{wildcards.burden_filter}/results/SKATO_{gene}\n")




checkpoint get_batches:
    params:
        nr_genes_batch = config["skato_batches"]
    input:
        burden_filter = "input_files/burden_filters/{burden_filter}.yml",
        gmat_list     = lambda wc: get_input_gene_output(wc)
    output:
        directory("run_folder/burden/{burden_filter}/batches")
    run:
        shell(f"mkdir -p run_folder/burden/{wildcards.burden_filter}/batches")
        prev_i = -1
        first = True
        with open(input.gmat_list, "r") as f:
            for i, g in enumerate(f.readlines()):
                if prev_i != int(i/params.nr_genes_batch):
                    prev_i = int(i/params.nr_genes_batch)
                    if not first:
                        w.close()
                    else:
                        first = False
                    w = open(f"run_folder/burden/{wildcards.burden_filter}/batches/{wildcards.burden_filter}_batch_{prev_i}", "w")
                w.write(g)


checkpoint get_variants_in_genes:
    """ This only works with genes as we assume they only exist on one chromosome, the other way would give memory problems """
    input:
        batch_list = "run_folder/burden/{burden_filter}/batches/{burden_filter}_batch_{n}",
        variants = lambda wc: expand("run_folder/burden/{{burden_filter}}/variant_matrix/vars/{gene}_var.tsv", gene=get_genes_in_batch(wc))
    output:
         directory("run_folder/burden/{burden_filter}/variant_matrix/annotated_vars/{n}")
    run:
        shell("mkdir -p run_folder/burden/{wildcards.burden_filter}/variant_matrix/annotated_vars/{wildcards.n}")
        for variants in input.variants:
            vars = False
            output_annotated = variants.split("/")
            output_annotated[4] = f"annotated_vars/{wildcards.n}"
            output_annotated[5] = output_annotated[5].replace(".tsv", "_distance.tsv")
            output_annotated = "/".join(output_annotated)
            with open(variants, "r") as f:
                lines = f.readlines()
                if len(lines) < 2:
                    shell(f"cp {variants} {output_annotated}")
                else:
                    vars = True
                    chrom = lines[1].split("\t")[1]
            if vars:
                distance_file = f"/proj/nobackup/sens2018106/analysis/pipeline_part3/input_files/exon_distance/chr{chrom}.npy"
                chrom_pos = np.load(distance_file, allow_pickle=True)
                chrom_pos = chrom_pos.item()
                if not chrom in chrom_pos:
                    chrom = "chr" + chrom
                max_len = len(chrom_pos[chrom]["pos"])
                with open(variants, "r") as f:
                    with open(output_annotated, "w") as w:
                        first = True
                        for line in f:
                            line = line.strip()
                            if first:
                                out_line = line + "\tDistance\tClosestElement\n"
                                first = False
                            else:
                                line_comp = line.split("\t")
                                pos = int(line_comp[2])
                                if max_len <= pos:
                                    dist = pos - max_len
                                    gene = chrom_pos[chrom]["gene"][-1]
                                else:
                                    dist = chrom_pos[chrom]["pos"][pos]
                                    gene = chrom_pos[chrom]["gene"][pos]
                                out_line = line + f"\t{dist}\t{gene}\n"
                            w.write(out_line)


checkpoint genotype_matrix:
    input:
        burden_filter = "input_files/burden_filters/{burden_filter}.yml",
        gene_tsv = lambda wc: get_gene_set(wc.burden_filter),
        vcf = config["filtered_vcf"]
    output:
        directory("run_folder/burden/{burden_filter}/variant_matrix")
    shell:
        """
        python3 scripts/get_genotype_matrix.py {input.vcf} {input.gene_tsv} {input.burden_filter}
        """
