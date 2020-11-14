# Author: Samvida S. Venkatesh
# Date: 24/06/20

PATH = [redacted]

library(TwoSampleMR)
library(tidyverse)

# Change these for each sensitivity analysis or as outcomes and 
# exposures change
exposures <- readRDS(paste(PATH, "/exposures/dfs_winner.rds", sep = ""))
outcomes <- readRDS(paste(PATH, "/outcomes/meta_outcomes.rds", sep = ""))
outcomes <- bind_rows(outcomes)

# Harmonise data ----

harmonised <- harmonise_data(exposures, outcomes)

# Perform MR ----

# Using Rucker's framework to choose method

# 1. IVW
IVW <- mr(harmonised, method_list = "mr_ivw")
IVW <- generate_odds_ratios(IVW)
# Adjust pval for multiple testing
IVW$pval <- p.adjust(IVW$pval, method = "fdr")
# Check heterogeneity
IVW_het <- mr_heterogeneity(harmonised, method_list = "mr_ivw")
# Merge 
IVW <- merge(IVW, IVW_het)

# 2. MR-Egger
egger <- mr(harmonised, method_list = "mr_egger_regression")
egger <- generate_odds_ratios(egger)
# Adjust pval for multiple testing
egger$pval <- p.adjust(egger$pval, method = "fdr")
# Check heterogeneity
egger_het <- mr_heterogeneity(harmonised, method_list = "mr_egger_regression")
# Merge 
egger <- merge(egger, egger_het)

# Merge with IVW results
all_MR <- IVW[, c("outcome", "exposure", "nsnp", "or", "or_lci95", "or_uci95", 
                  "pval", "Q", "Q_pval")]
colnames(all_MR)[4:9] <- c("OR_IVW", "LCI_IVW", "UCI_IVW", "pval_IVW", 
                           "Q_IVW", "Qpval_IVW")
all_MR <- merge(all_MR, egger[, c("outcome", "exposure", "nsnp", "or", "or_lci95",
                                  "or_uci95",  "pval", "Q")])
colnames(all_MR)[10:14] <- c("OR_Egger", "LCI_Egger", "UCI_Egger", 
                             "pval_Egger", "Q_Egger")

# Check if MR-Egger is significantly less heterogeneous

all_MR$diff_Q <- all_MR$Q_IVW - all_MR$Q_Egger
all_MR$Qpval_egger <- pchisq(all_MR$diff_Q, df = 1)

# Method to use
all_MR$method <- ifelse(all_MR$Qpval_IVW >= 0.05, "IVW", "Egger")

# 3. Weighted median for sensitivity
wtmed <- mr(harmonised, method_list = "mr_weighted_median")
wtmed <- generate_odds_ratios(wtmed)
# Adjust pval for multiple testing
wtmed$pval <- p.adjust(wtmed$pval, method = "fdr")

all_MR <- merge(all_MR, wtmed[, c("outcome", "exposure", "or", "or_lci95", "or_uci95",
                                  "pval")])
colnames(all_MR)[18:21] <- c("OR_weighted_median", "LCI_weighted_median", 
                             "UCI_weighted_median", "pval_weighted_median")

write.table(all_MR, paste(PATH, 
                          "/results/obesity_meta_wrh_mr_results_280720.txt", 
                          sep = ""),
            sep = "\t",
            quote = F, row.names = F)