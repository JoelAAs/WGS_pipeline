# Script modified from http://cnsgenomics.com/software/gcta/#GREMLinWGSorimputeddata
args <- commandArgs(trailingOnly = T)

ld_score_file = args[1]
association   = args[2]
path = paste("run_folder", "GCTA", "tests", association, sep="/")

lds_seg = read.table(ld_score_file, header=T,colClasses=c("character",rep("numeric",8)))
quartiles=summary(lds_seg$ldscore_SNP)

lb1 = which(lds_seg$ldscore_SNP <= quartiles[2]) 
lb2 = which(lds_seg$ldscore_SNP > quartiles[2] & lds_seg$ldscore_SNP <= quartiles[3])
lb3 = which(lds_seg$ldscore_SNP > quartiles[3] & lds_seg$ldscore_SNP <= quartiles[5])
lb4 = which(lds_seg$ldscore_SNP > quartiles[5])

multi_grm_file <- file(
  paste0(path, "/", association, "_multi_grm.txt"),
  "w+"
)
group_idx = c(lb1, lb2, lb3, lb4)
group_id = 0
for (snp_idx in group_idx) {
  group_id = group_id + 1
  name = paste0("snp_group", group_id)
  write.table(lds_seg$SNP[snp_idx], 
              paste0(path, "/", name, ".txt"),
              row.names=F, quote=F, col.names=F)
  write(name, multi_grm_file)
}
close(multi_grm_file)
