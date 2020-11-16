# Author: Samvida S. Venkatesh
# Date: 04/08/20

PATH = [redacted]

library(TwoSampleMR)
library(dplyr)

instruments <- readRDS(paste(PATH, 
                             "/exposures/separate_MVMR_instruments.rds", 
                             sep = ""))
outcome_names <- c("endometriosis", "HMB", "pcos", "pre_or_eclamps", "uf")

outcomes <- lapply(instruments, function (ins) {
  # Go through each outcome to get relevant SNPs for each set of MVMR instruments
  res <- lapply(outcome_names, function (oname) {
    # Read data in 2SMR format
    tmp <- read_outcome_data(filename = paste(PATH, "/meta_outcomes/meta_", 
                                       oname, ".txt", sep = ""),
                      snps = ins$SNP, 
                      sep = "\t", 
                      snp_col = "MarkerName", beta_col = "beta", se_col = "se",
                      effect_allele_col = "Allele1", other_allele_col = "Allele2",
                      eaf_col = "Freq1", pval_col = "P.value")
    tmp$outcome <- oname
    return (tmp)
  })
  res <- bind_rows(res)
  return (res)
})

saveRDS(outcomes, paste(PATH, "/outcomes/MVMR_outcomes.rds", sep = ""))