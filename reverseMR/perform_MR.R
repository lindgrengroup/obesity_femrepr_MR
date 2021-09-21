# Author: Samvida S. Venkatesh
# Date: 06/04/20

library(TwoSampleMR)
library(tidyverse)

# Construct reproductive disease instruments ----

repr_traits <- c("endometriosis", "pcos", "uf")

## Function to flip SNPs whose effect sizes are negative
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

## Read instrument data 
repr_instr <- lapply(repr_traits, function (r) {
  res <- read_exposure_data(filename = paste0("instruments/", r, "_formatted.txt"),
                            sep = "\t",
                            phenotype_col = "DISEASE", snp_col = "SNP",
                            beta_col = "BETA", se_col = "S.E.",
                            effect_allele_col = "RISK ALLELE", 
                            other_allele_col = "OTHER ALLELE", 
                            eaf_col = "RISK ALLELE FREQUENCY",
                            pval_col = "P-VALUE",
                            ncase_col = "NCASES", ncontrol_col = "NCONTROLS")
  res <- flipSNPs(res)
  return (res)
})
names(repr_instr) <- repr_traits

## Assess instrument strength (F-statistics)

calcFStatistic <- function(x) {
  # formula for F-statistic: 1/L * sum over all SNPs(beta.exp^2 / var.exp)
  num <- (x$beta.exposure)^2
  denom <- (x$se.exposure)^2
  fsum <- sum(num/denom, na.rm = T)
  fmean <- fsum / length(num)
  return (fmean)
}

# Calculate dilution percentage for each constructed instrument
lapply(repr_instr, calcFStatistic)

repr_instr <- bind_rows(repr_instr)

saveRDS(repr_instr, "instruments/repr_trait_instruments.rds")

# Get obesity outcome data ----

obesity_traits <- c("bmi", "whr", "whradjbmi")
instrument_SNPs <- unique(repr_instr$SNP)

outcome_data <- lapply(obesity_traits, function (o) {
  df <- read_outcome_data(
    snps = instrument_SNPs,
    filename = paste0("/well/lindgren/resources/GIANT_UKB_2018_BMI_WHR_WHRadjBMI/",
                      o, ".giant-ukbb.meta-analysis.females.23May2018.txt.gz"),
    sep = " ",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "Tested_Allele",
    other_allele_col = "Other_Allele",
    eaf_col = "Freq_Tested_Allele",
    pval_col = "P",
    samplesize_col = "N"
  )
  df$outcome <- o
  return (df)
})
names(outcome_data) <- obesity_traits
outcome_data <- bind_rows(outcome_data)

saveRDS(outcome_data, "outcomes/obesity_outcomes.rds")

# Harmonise data ----

h_data <- harmonise_data(repr_instr, obesity_outcomes)

# Perform MR ----

exposures <- repr_instr
outcomes <- obesity_outcomes
harmonised <- h_data

# Using Rucker's framework to choose method

# 1. IVW
IVW <- mr(harmonised, method_list = "mr_ivw")
# Adjust pval for multiple testing
# IVW$pval <- p.adjust(IVW$pval, method = "fdr")
# Check heterogeneity
IVW_het <- mr_heterogeneity(harmonised, method_list = "mr_ivw")
# Merge 
IVW <- merge(IVW, IVW_het)

# 2. MR-Egger
egger <- mr(harmonised, method_list = "mr_egger_regression")
# Adjust pval for multiple testing
# egger$pval <- p.adjust(egger$pval, method = "fdr")
# Check heterogeneity
egger_het <- mr_heterogeneity(harmonised, method_list = "mr_egger_regression")
# Merge 
egger <- merge(egger, egger_het)

# Merge with IVW results
all_MR <- IVW[, c("outcome", "exposure", "nsnp", "b", "se", "pval", 
                  "Q", "Q_pval")]
colnames(all_MR)[4:8] <- c("beta_IVW", "SE_IVW", "pval_IVW", 
                           "Q_IVW", "Qpval_IVW")
all_MR <- merge(all_MR, egger[, c("outcome", "exposure", "nsnp", "b", "se",
                                  "pval", "Q")])
colnames(all_MR)[9:12] <- c("beta_Egger", "SE_Egger", "pval_Egger", "Q_Egger")

# Check if MR-Egger is significantly less heterogeneous

all_MR$diff_Q <- all_MR$Q_IVW - all_MR$Q_Egger
all_MR$Qpval_egger <- pchisq(all_MR$diff_Q, df = 1)

# Method to use
all_MR$method <- ifelse(all_MR$Qpval_IVW >= 0.05, "IVW", 
                        ifelse(all_MR$Qpval_egger >= 0.05, "Egger", "WTMed"))

# 3. Weighted median for sensitivity
wtmed <- mr(harmonised, method_list = "mr_weighted_median")
# Adjust pval for multiple testing
# wtmed$pval <- p.adjust(wtmed$pval, method = "fdr")

all_MR <- merge(all_MR, wtmed[, c("outcome", "exposure", "b", "se", "pval")])
colnames(all_MR)[16:18] <- c("beta_weighted_median", "SE_weighted_median", 
                             "pval_weighted_median")

write.table(all_MR, "results/mr_results_120421.txt",
            sep = "\t",
            quote = F, row.names = F)

# Steiger filtering of SNPs ----

# Calculate r.exposure assuming ~10% prevalence of all 
# reproductive disorders, i.e. PCOS, UF, and endometriosis
harmonised$r.exposure <- get_r_from_lor(lor = harmonised$beta.exposure,
                            af = harmonised$eaf.exposure,
                            ncase = harmonised$ncase.exposure,
                            ncontrol = harmonised$ncontrol.exposure,
                            prevalence = 0.1)
harmonised$r.outcome <- get_r_from_bsen(b = harmonised$beta.outcome,
                                             se = harmonised$se.outcome,
                                             n = harmonised$samplesize.outcome)

# Steiger filtering
steiger_res <- steiger_filtering(harmonised)