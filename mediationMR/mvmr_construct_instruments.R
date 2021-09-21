# Author: Samvida S. Venkatesh
# Date: 04/08/20

# List instrument SNPs for each exposure and mediator ----
# Hormone traits
hormones <- readRDS("mediators/mediators_as_instruments.rds")
instrument_SNPs <- hormones[, c("SNP", "exposure")]

# Obesity traits
SNPlist <- lapply(c("bmi", "whr", "whradjbmi"), function (x) {
  snps <- read.table(paste("exposures/", x, "_indexsnps_females.txt", sep = ""), 
                     header = T, sep = " ", stringsAsFactors = F)$SNP
  snps <- sub("\\:.*", "", snps)
  return (data.frame(SNP = snps, exposure = x))
})

instrument_SNPs <- bind_rows(instrument_SNPs, bind_rows(SNPlist))

saveRDS(instrument_SNPs, "exposures/all_instrument_SNPs.rds")

# Construct MVMR instruments for each combination of obesity x hormone ----

# First collect all SNPs for hormone traits

hormones <- c("Leptin", "Fasting insulin", "Insulin sensitivity",
              "Age at menarche", "Age at menopause")
hormone_codes <- c("ebi-a-GCST003367", "ebi-a-GCST005185", "ebi-a-GCST005178",
                   "ieu-a-1095", "ieu-a-1004")
names(hormone_codes) <- hormones

hormone_MVMR <- lapply(hormones, function (h) {
  # Get values for given SNPs
  res <- extract_outcome_data(snps = instrument_SNPs$SNP, 
                              outcomes = hormone_codes[h],
                              proxies = F)
  # Convert to exposure data format
  common_cols <- c("id", "effect_allele", "other_allele", "eaf",
                   "beta", "se", "pval")
  cols_to_keep <- c("SNP", "outcome", paste(common_cols, "outcome", sep = "."))
  cols_rename <- c("SNP", "exposure", paste(common_cols, "exposure", sep = "."))
  
  res <- res[, cols_to_keep]
  colnames(res) <- cols_rename
  res$exposure <- h
  
  return (res)
  
})

saveRDS(hormone_MVMR, "mediators/MVMR_instruments.rds")

# Now collect all SNPs for obesity traits

obesity <- c("bmi", "whr", "whradjbmi")

obesity_MVMR <- lapply(obesity, function (o) {
  # Get values for given SNPs
  res <- read_outcome_data(filename = paste("/well/lindgren/samvida/PULIT_2018/fat-distn.giant.ukbb.meta-analysis.", 
                                            o, ".females.tbl.txt", sep = ""),
                           snps = instrument_SNPs$SNP, 
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

MVMR_instruments <- bind_rows(bind_rows(hormone_MVMR), bind_rows(obesity_MVMR))

saveRDS(MVMR_instruments, "exposures/MVMR_instruments.rds")

# Subset instruments for obesity x hormone combinations of traits ----

instrument_SNPs <- readRDS("exposures/all_instrument_SNPs.rds")
MVMR_instruments <- readRDS("exposures/MVMR_instruments.rds")

hormones <- c("Leptin", "Fasting insulin", "Insulin sensitivity",
              "Age at menarche", "Age at menopause")
obesity <- c("bmi", "whr", "whradjbmi")

split_MVMR_instruments <- list()

for (h in hormones) {
  for (o in obesity) {
    # Get relevant instrument SNPs
    SNPs <- instrument_SNPs[instrument_SNPs$exposure %in% c(h, o), "SNP"]
    # Get instrument details for those SNPs
    subset_MVMR <- subset(MVMR_instruments, MVMR_instruments$SNP %in% SNPs &
                            MVMR_instruments$exposure %in% c(h, o))
    split_MVMR_instruments[[paste(o, h, sep = "_")]] <- subset_MVMR
  }
}

saveRDS(split_MVMR_instruments, "exposures/separate_MVMR_instruments.rds")

# Subset instruments for obesity x 2 hormone combinations ----

instrument_SNPs <- readRDS("exposures/all_instrument_SNPs.rds")
MVMR_instruments <- readRDS("exposures/MVMR_instruments.rds")

obesity <- c("bmi", "whr", "whradjbmi")

joint_leptin_FI_instruments <- list()
joint_leptin_ISI_instruments <- list()

for (o in obesity) {
  # Leptin and FI SNPs with obesity trait
  lep_FI_SNPs <- 
    instrument_SNPs[instrument_SNPs$exposure %in% 
                      c("Leptin", "Fasting insulin", o), 
                    "SNP"]
  # Leptin and ISI SNPs with obesity trait
  lep_ISI_SNPs <- 
    instrument_SNPs[instrument_SNPs$exposure %in% 
                      c("Leptin", "Insulin sensitivity", o), 
                    "SNP"]
  
  # Get instrument details for those SNPs
  subset_lep_FI <- subset(MVMR_instruments, 
                          MVMR_instruments$SNP %in% lep_FI_SNPs &
                            MVMR_instruments$exposure 
                          %in% c("Leptin", "Fasting insulin", o))
  joint_leptin_FI_instruments[[o]] <- subset_lep_FI
  
  subset_lep_ISI <- subset(MVMR_instruments, 
                          MVMR_instruments$SNP %in% lep_ISI_SNPs &
                            MVMR_instruments$exposure 
                          %in% c("Leptin", "Insulin sensitivity", o))
  joint_leptin_ISI_instruments[[o]] <- subset_lep_ISI
}

saveRDS(joint_leptin_FI_instruments, "exposures/joint_lep_FI_instruments.rds")
saveRDS(joint_leptin_ISI_instruments, "exposures/joint_lep_ISI_instruments.rds")


