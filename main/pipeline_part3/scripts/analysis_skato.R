.libPaths("/proj/sens2018106/softwares/R_packages")

library(SKAT)
library(yaml)

args <- commandArgs(trailingOnly = T)
#setwd(dirname(args[1]))

## Arguments

variant_ann_filename      <- args[1]
phenotype_vector_filename <- args[2]
config_file               <- args[3]
chosen_tests_file         <- args[4]
excluded_samples          <- args[5]
covar_files               <- args[6:length(args)]


genotype_matrix_filename   <- gsub("var_distance.tsv", "gmat.tsv", variant_ann_filename, fixed = T)
genotype_matrix_filename   <-  gsub("annotated_vars/[0-9]+", "gmats", genotype_matrix_filename, perl = T)
variant_ann_filename_c <- gsub("_var_distance.tsv", "_var_final.tsv", variant_ann_filename, fixed = T)
geneset_name           <- gsub("_gmat.tsv", "", tail(strsplit(genotype_matrix_filename, "/")[[1]], 1), fixed = T)
output_filename        <- paste0(gsub("variant_matrix/gmats","results",dirname(genotype_matrix_filename),),
                            paste("/SKATO", geneset_name, "result.tsv", sep="_"))


## Load genotype matrix and variant annotation
genotype_matrix  <- read.csv(genotype_matrix_filename, sep = "\t", stringsAsFactors = F, check.names = F)
variant_ann      <- read.csv(variant_ann_filename, sep = "\t", stringsAsFactors = F, check.names = F)
phenotype        <- read.csv(phenotype_vector_filename, sep ="\t", stringsAsFactors = F, check.names = F)
excluded_df <- read.csv(excluded_samples, stringsAsFactors = F, check.names = F, header = F)
filter_config    <- read_yaml(config_file)
row.names(phenotype) <- phenotype$IID


get_hwe <- function(curr_phenotype_group) {
  df <- read.csv(
    paste0("run_folder/hwe/", curr_phenotype_group,"_snps"),
    sep = "\t", stringsAsFactors = F, header = F)
  return(df)
}

get_vf_for_group <- function(phenotype_groups, phenotype, variant_ann, genotype_matrix) {
  for (curr_phenotype_group in phenotype_groups) {
    samples <- row.names(phenotype[phenotype[, curr_phenotype_group] != -9, ])
    samples <- samples[samples %in% row.names(genotype_matrix)]

    variant_ann[, paste0(curr_phenotype_group, "_COUNT")] <- sapply(
      variant_ann$PK, function (x) sum(genotype_matrix[samples, as.character(x)], na.rm = T))

    variant_ann[, paste0(curr_phenotype_group, "_VF")] <- sapply(
      variant_ann[, paste0(curr_phenotype_group, "_COUNT")], function(x) x/(length(samples)*2))  # variant freq
   }
  return(variant_ann)
}

read_covar <- function(covar_file_paths, current_assoc, samples){
  covar_file    <- covar_file_paths[grepl(paste0(".*", current_assoc,".*"), covar_file_paths)][1]
  covar_csv     <- read.csv(covar_file, sep="\t", stringsAsFactors = F, check.names = F)
  covar_csv$IID <- covar_csv[, 1]
  row.names(covar_csv) <- covar_csv$IID
  covar_matrix <- as.matrix(covar_csv[samples, 2:4])
  return(covar_matrix)
}
## Get samples
geno_sample_names         <- colnames(genotype_matrix)[-1] #Remove PK
geno_sample_names         <- geno_sample_names[!(geno_sample_names %in% excluded_df$V1)]
variant_ann$COUNT         <- sapply(variant_ann$PK, function (x)
  sum(genotype_matrix[genotype_matrix$PK == x, geno_sample_names], na.rm = T))
variant_ann$VF            <- sapply(variant_ann$COUNT, function(x) x/(length(geno_sample_names)*2))  # variant freq

## Get weights
sg <- function(x, h, mid) {
  1 - 1/(1+exp(h*(mid-x)))
  }



variant_ann$Distance <- as.numeric(variant_ann$Distance)
variant_ann$Weight   <- sapply(variant_ann$Distance, function(x)
  round(sg(x, 0.2, 40), digits=3))
variant_ann$hweid <- paste(paste(variant_ann$CHROM, variant_ann$POS, sep="_"), variant_ann$ALT, variant_ann$REF, sep=":")


## Run filters
prev_bool_idx <- rep(T, nrow(variant_ann))

for (filter_type in c("all", "not", "is", "less", "greater")) {
  current_filter_set <- filter_config$burden_filters[[filter_type]]
  filter_targets <- names(current_filter_set)
  for (target in filter_targets) {
    values <- current_filter_set[[target]]
    switch(filter_type,
           "all" = rep(T, nrow(variant_ann)),
           "not" = !(variant_ann[, target] %in% values),
           "is" = variant_ann[, target] %in% values,
           "less" = as.numeric(variant_ann[, target]) < values,
           "greater" = as.numeric(variant_ann[, target]) > values,
           stop(paste0("Unkown filtertype: '", filter_type, "'"))
    ) -> bool_index
    bool_index[is.na(bool_index)] <- F
    prev_bool_idx  <- prev_bool_idx & bool_index
  }
}
prev_bool_idx <- prev_bool_idx & variant_ann$COUNT > 1  # remove singletons
if (nrow(variant_ann) != 0){
  prev_bool_idx <- prev_bool_idx & variant_ann$Weight # remove 0 weighted variants
}


## Can we assume that the Genotype matrix and variant annotation allways have the same order? They should!
print(paste0("Keeping ", sum(prev_bool_idx), " variants out of: ", length(prev_bool_idx)))
if (sum(prev_bool_idx) < 2) {
  out_file <- file(output_filename, "w+")
  close(out_file)
  # R saknar Stop(exit_success) så det blir if-else ist
} else {

  weights                   <- variant_ann$Weight[prev_bool_idx]
  phenotype                 <- phenotype[phenotype$IID %in% geno_sample_names, ]

  genotype_matrix           <- genotype_matrix[prev_bool_idx, ]
  rownames(genotype_matrix) <- genotype_matrix$PK


  genotype_matrix$PK        <- NULL
  if (!all(geno_sample_names == phenotype$IID)){
    warning("The genotypematrix sample order is different from the phenotype data.. reordering")
    genotype_matrix = genotype_matrix[, phenotype$IID]
  }

  write.table(variant_ann, file = gsub(".tsv", "unfiltered.tsv", variant_ann_filename, fixed = T), sep="\t")
  genotype_matrix <- t(as.matrix(genotype_matrix)) # Make matrix and transpose
  variant_ann <- variant_ann[variant_ann$PK %in% colnames(genotype_matrix), ]
  phenotype_groups <- readLines(chosen_tests_file)
  variant_ann <- get_vf_for_group(phenotype_groups, phenotype, variant_ann, genotype_matrix)
  write.csv2(variant_ann, variant_ann_filename_c)


  out_file <- file(output_filename, "w+")
  for (curr_phenotype_group in phenotype_groups){
    print(curr_phenotype_group)
    n_cases <- length(row.names(phenotype[phenotype[, curr_phenotype_group] == 2, ]))
    threshold_VF <- (1/n_cases)^(1/2)

    vf_pks  <- variant_ann[variant_ann[, paste0(curr_phenotype_group, "_VF")] > threshold_VF, "PK"]

    hwe_df   <- get_hwe(curr_phenotype_group)
    hwe_pks  <- variant_ann[variant_ann$hweid %in% hwe_df$V2, "PK"]
    pass_pks <- colnames(genotype_matrix)
    pass_pks <- pass_pks[!pass_pks %in% hwe_pks]
    pass_pks <- pass_pks[!pass_pks %in% vf_pks]
    samples <- row.names(phenotype[phenotype[, curr_phenotype_group] != -9, ])
    samples <- samples[samples %in% row.names(genotype_matrix)]
    if (length(samples) > 500) {
      covar_matrix <- read_covar(covar_files, curr_phenotype_group, samples)
      obj <- SKAT_Null_Model((phenotype[samples, curr_phenotype_group]-1) ~ covar_matrix, out_type="D")
      results <- SKAT(
        genotype_matrix[samples, pass_pks],
        obj,
        impute.method = "bestguess",
        method = "SKATO",
        kernel = "linear.weighted",
        weights.beta = weights
      )
      write(paste(geneset_name, curr_phenotype_group, results$p.value, results$param$n.marker.test, sep="\t"), out_file)
    }
  }
  close(out_file)
}
