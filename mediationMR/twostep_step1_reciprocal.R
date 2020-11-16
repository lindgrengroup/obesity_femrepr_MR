# Author: Samvida S. Venkatesh
# Date: 03/08/2020

PATH = [redacted]

library(TwoSampleMR)
library(dplyr)

# Build instruments ----

leptin <- extract_instruments(outcomes = "ebi-a-GCST003367", p1 = 5e-06,
                              clump = T)
leptin$exposure <- "Leptin"

fasting_insulin <- extract_instruments(outcomes = "ebi-a-GCST005185", p1 = 5e-06,
                                       clump = T)
fasting_insulin$exposure <- "Fasting insulin"

insulin_sensitivity <- extract_instruments(outcomes = "ebi-a-GCST005178", p1 = 5e-06,
                                           clump = T)
insulin_sensitivity$exposure <- "Insulin sensitivity"

mediator_instruments <- list(leptin, insulin_sensitivity, fasting_insulin)

# Flip SNPs 
# SNPs with negative effect sizes need to be flipped so they have a 
# positive effect size; this involves:
# 1. reversing the sign of beta
# 2. swapping effect allele and other allele
# 3. changing effect allele frequency to 1 - eaf

flipSNPs <- function(x) {
  # get SNPs whose effect sizes are negative
  indexes <- which(x$beta.exposure < 0)
  # reverse the sign of beta
  x[indexes, "beta.exposure"] <- -x[indexes, "beta.exposure"]
  # swap effect allele and other allele
  orig_A1 <- x[indexes, "effect_allele.exposure"]
  orig_A2 <- x[indexes, "other_allele.exposure"]
  x[indexes, "effect_allele.exposure"] <- orig_A2
  x[indexes, "other_allele.exposure"] <- orig_A1
  # change eaf to 1 - eaf
  x[indexes, "eaf.exposure"] <- 1 - x[indexes, "eaf.exposure"]
  
  return (x)
}

mediator_instruments <- lapply(mediator_instruments, flipSNPs)

# Assess instrument strength (F-statistics) ----

calcFStatistic <- function(x) {
  # formula for F-statistic: 1/L * sum over all SNPs(beta.exp^2 / var.exp)
  num <- (x$beta.exposure)^2
  denom <- (x$se.exposure)^2
  fsum <- sum(num/denom, na.rm = T)
  fmean <- fsum / length(num)
  return (fmean)
}

# Calculate dilution percentage for each constructed instrument
lapply(mediator_instruments, calcFStatistic)

mediator_instruments <- bind_rows(mediator_instruments)
saveRDS(mediator_instruments, 
        paste(PATH, "/mediators/mediators_as_instruments.rds", sep = ""))

# Collect exposure data (exposures as outcomes) ----

all_SNPs <- unique(mediator_instruments$SNP)

# Modify exposure files SNP column from rsid:allele1:allele2 to rsid
exposure_files <- list.files("/well/lindgren/samvida/PULIT_2018", 
                             pattern = "*.tbl")

lapply(exposure_files, function (fname) {
  tmp <- read.table(paste("/well/lindgren/samvida/PULIT_2018/", 
                          fname, sep = ""), 
                    sep = "\t", header = T, comment.char = "¬", 
                    stringsAsFactors = F)
  tmp$SNP <- sub("\\:.*", "", tmp$SNP)
  write.table(tmp, paste("/well/lindgren/samvida/PULIT_2018/", 
                         fname, ".txt", sep = ""), 
              sep = "\t", quote = F, row.names = F)
})

exposure_files <- list.files("/well/lindgren/samvida/PULIT_2018", 
                             pattern = "*.txt")
ename <- c("bmi", "whr", "whradjbmi")

exposures <- lapply(1:3, function (i) {
  # Read data in 2SMR format
  res <- read_outcome_data(filename = 
                             paste("/well/lindgren/samvida/PULIT_2018/", 
                                   exposure_files[i], sep = ""),
                           snps = all_SNPs, 
                           sep = "\t", 
                           snp_col = "SNP", beta_col = "Effect", 
                           se_col = "StdErr",
                           effect_allele_col = "Allele1", 
                           other_allele_col = "Allele2",
                           eaf_col = "Freq1", pval_col = "P-value")
  res$outcome <- ename[i]
})

exposures <- bind_rows(exposures)
saveRDS(exposures, paste(PATH, "/exposures/exposures_as_outcomes.rds", 
                         sep = ""))

# Harmonise ----

set1 <- readRDS(paste(PATH, "/mediators/mediators_as_instruments.rds", 
                      sep = ""))
set2 <- readRDS(paste(PATH, "/exposures/exposures_as_outcomes.rds", 
                      sep = ""))

harmonised <- harmonise_data(set1, set2)

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

write.table(all_MR, paste(PATH, 
                          "/results/exposure_mediator_reciprocal_030820.txt", 
                          sep = ""), sep = "\t", quote = F, row.names = F)

# Directionality: MR-Steiger ----

harmonised$samplesize.outcome <- ifelse(harmonised$outcome == "bmi", 434794, 
                                        ifelse(harmonised$outcome == "whr", 381152,
                                               379501))
dir <- directionality_test(harmonised)

