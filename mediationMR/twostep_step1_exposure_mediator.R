# Author: Samvida S. Venkatesh
# Date: 03/08/2020

PATH = [redacted]

library(TwoSampleMR)
library(dplyr)

# Exposures ----

# Get exposures from previously calculated instruments for BMI, WHR, and WHRadjBMI
# Female-specific SNPs with female-specific weights
exposures <- readRDS(paste(PATH, "/exposures/dfs_winner.rds", sep = ""))

# Get mediator outcomes ----

leptin <- extract_outcome_data(snps = exposures$SNP,
                               outcomes = "ebi-a-GCST003367")

fasting_insulin <- extract_outcome_data(snps = exposures$SNP,
                           outcomes = "ebi-a-GCST005185")

insulin_sensitivity <- extract_outcome_data(snps = exposures$SNP,
                                            outcomes = "ebi-a-GCST005178")

mediators <- bind_rows(leptin, fasting_insulin, insulin_sensitivity)

saveRDS(mediators, paste(PATH, "/mediators/mediator_outcomes.rds", sep = ""))

# Harmonise data ----

harmonised <- harmonise_data(exposures, mediators)

# Perform MR by picking best method as outlined by Rucker's framework ----

# 1. IVW
IVW <- mr(harmonised, method_list = "mr_ivw")
# Adjust pval for multiple testing
IVW$pval <- p.adjust(IVW$pval, method = "fdr")
# Check heterogeneity
IVW_het <- mr_heterogeneity(harmonised, method_list = "mr_ivw")
# Merge 
IVW <- merge(IVW, IVW_het)

# 2. MR-Egger
egger <- mr(harmonised, method_list = "mr_egger_regression")
# Adjust pval for multiple testing
egger$pval <- p.adjust(egger$pval, method = "fdr")
# Check heterogeneity
egger_het <- mr_heterogeneity(harmonised, method_list = "mr_egger_regression")
# Merge 
egger <- merge(egger, egger_het)

# Merge with IVW results
all_MR <- IVW[, c("outcome", "exposure", "nsnp", "b", "se", 
                  "pval", "Q", "Q_pval")]
colnames(all_MR)[4:8] <- c("beta_IVW", "se_IVW", "pval_IVW", "Q_IVW", "Qpval_IVW")
all_MR <- merge(all_MR, egger[, c("outcome", "exposure", "nsnp", "b", "se",
                                  "pval", "Q")])
colnames(all_MR)[9:12] <- c("beta_Egger", "se_Egger", "pval_Egger", "Q_Egger")

# Check if MR-Egger is significantly less heterogeneous

all_MR$diff_Q <- all_MR$Q_IVW - all_MR$Q_Egger
all_MR$Qpval_egger <- pchisq(all_MR$diff_Q, df = 1)

# Method to use
all_MR$method <- ifelse(all_MR$Qpval_IVW >= 0.05, "IVW", 
                        ifelse(all_MR$Qpval_egger >= 0.5, "Egger", "WTMed"))

# 3. Weighted median for sensitivity
wtmed <- mr(harmonised, method_list = "mr_weighted_median")
# Adjust pval for multiple testing
wtmed$pval <- p.adjust(wtmed$pval, method = "fdr")

all_MR <- merge(all_MR, wtmed[, c("outcome", "exposure", "b", "se", "pval")])
colnames(all_MR)[16:18] <- c("beta_weighted_median", "se_weighted_median", 
                             "pval_weighted_median")

write.table(all_MR, paste(PATH, "/results/exposure_mediator_030820.txt", 
                          sep = ""), sep = "\t", quote = F, row.names = F)

# Directionality: MR-Steiger ----

dir <- directionality_test(harmonised)