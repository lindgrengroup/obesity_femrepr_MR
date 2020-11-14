# Author: Samvida S. Venkatesh
# Date: 24/06/20

library(TwoSampleMR, lib.loc = "/well/lindgren/samvida/R/lib")

outcome_files <- list.files("/well/lindgren/samvida/FINNGEN_R3")

outcome_names <- c("Endometriosis", "Excessive or irregular menstruation", 
                   "Infertility", "Miscarriage", "PCOS", "Preeclampsia or eclampsia",
                   "Uterine fibroids")
  
# Write data to format that can be read by 2SMR
lapply(1:length(outcome_names), function (i) {
  
  res <- read.table(paste("/well/lindgren/samvida/FINNGEN_R3/", outcome_files[i], 
                          sep = ""),
                    sep = "\t", header = T, comment.char = "?")
  res$phenotype <- outcome_names[i]
  
  write.table(res, paste("/well/lindgren/samvida/obesity_wrh/twosampleMR/finngen_outcomes/", 
                         outcome_files[i], ".txt", sep = ""), 
              sep = "\t", row.names = F, quote = F)
})

# Extract outcome information for all SNPs present in all instruments; 
# then subset the dataset as needed for each MR

exposures <- readRDS("/well/lindgren/samvida/obesity_wrh/twosampleMR/exposures/dfs_winner.rds")
# Get list of SNPs to extract information for
all_SNPs <- unique(exposures$SNP)

## Read FinnGen data ----

outcome_files <- list.files("/well/lindgren/samvida/obesity_wrh/twosampleMR/finngen_outcomes", 
                            pattern = "finngen_r3_*")
  
outcomes <- lapply(outcome_files, function (fname) {
  # Read data in 2SMR format
  read_outcome_data(filename = paste("/well/lindgren/samvida/obesity_wrh/twosampleMR/finngen_outcomes/", 
                                     fname, sep = ""),
                    snps = all_SNPs, 
                    sep = "\t", 
                    phenotype_col = "phenotype", 
                    snp_col = "rsids", beta_col = "beta", se_col = "sebeta",
                    effect_allele_col = "alt", other_allele_col = "ref",
                    eaf_col = "maf", pval_col = "pval", gene_col = "nearest_genes")
})

saveRDS(outcomes, "/well/lindgren/samvida/obesity_wrh/twosampleMR/outcomes/finngen_r3_outcomes.rds")
