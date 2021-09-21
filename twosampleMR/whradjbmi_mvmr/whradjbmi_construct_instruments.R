# Author: Samvida S. Venkatesh
# Date: 21/09/21

library(tidyverse)
library(TwoSampleMR)

# Get instrument SNPs for both WHR and BMI ----

instr <- readRDS("exposures/dfs_winner.rds")
whr_and_bmi_SNPs <- unique(subset(instr, 
                                  instr$exposure != "whradjbmi")$SNP)

# Get effect estimates ----

for_mvmr <- lapply(c("bmi", "whr"), function (o) {
  # Get values for given SNPs
  res <- read_outcome_data(filename = paste("/well/lindgren/samvida/PULIT_2018/fat-distn.giant.ukbb.meta-analysis.", 
                                            o, ".females.tbl.txt", sep = ""),
                           snps = whr_and_bmi_SNPs, 
                           sep = "\t", 
                           snp_col = "SNP", beta_col = "Effect", se_col = "StdErr",
                           effect_allele_col = "Allele1", other_allele_col = "Allele2",
                           eaf_col = "Freq1", pval_col = "P-value")
  # Convert to exposure data format
  common_cols <- c("id", "effect_allele", "other_allele", "eaf",
                   "beta", "se", "pval")
  cols_to_keep <- c("SNP", "outcome", paste(common_cols, "outcome", sep = "."))
  cols_rename <- c("SNP", "exposure", paste(common_cols, "exposure", sep = "."))
  
  res <- res[, cols_to_keep]
  colnames(res) <- cols_rename
  res$exposure <- o
  
  return (res)
})

for_mvmr <- bind_rows(bind_rows(for_mvmr))

saveRDS(for_mvmr, "exposures/mvmr_whr_and_bmi_instr.rds")

