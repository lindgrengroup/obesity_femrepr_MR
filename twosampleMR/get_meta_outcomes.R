# Author: Samvida S. Venkatesh
# Date: 02/07/20

PATH = [redacted]

library(TwoSampleMR)

# Extract outcome information for all SNPs present in all instruments; 
# then subset the dataset as needed for each MR

exposures <- readRDS(paste(PATH, "/exposures/dfs_winner.rds", sep = ""))
# Get list of SNPs to extract information for
all_SNPs <- unique(exposures$SNP)

## Read FinnGen data ----

outcome_files <- list.files(paste(PATH, "/meta_outcomes", sep = ""), 
                            pattern = "*.txt")
outcome_names <- c("Endometriosis", "Excessive or frequent menstruation",
                   "Infertility", "Miscarriage", "PCOS", 
                   "Preeclampsia or eclampsia", "Uterine fibroids")

outcomes <- lapply(1:length(outcome_files), function (i) {
  # Read data in 2SMR format
  df <- read_outcome_data(filename = paste(PATH, "/meta_outcomes/", 
                                           outcome_files[i], sep = ""),
                          snps = all_SNPs, 
                          sep = "\t", 
                          snp_col = "MarkerName", beta_col = "beta", se_col = "se",
                          effect_allele_col = "Allele1", other_allele_col = "Allele2",
                          eaf_col = "Freq1", pval_col = "P.value")
  df$outcome <- outcome_names[i]
  return (df)
})

saveRDS(outcomes, paste(PATH, "/outcomes/meta_outcomes.rds", sep = ""))
