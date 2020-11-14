# Author: Samvida S. Venkatesh
# Date: 16/01/2020

# Load relevant packages
library(tidyverse)     # for data manipulation
library(TwoSampleMR)   # to perform two-sample MR
library(MRInstruments) # to get GWAS Catalog data for MR

# Construct Instruments ----
### Read exposure data for WHR, BMI, WHRadjBMI ----

# Read Pulit et al. data
# Name all the data by trait and whether the SNPs are from combined, male, or 
# female analyses

traits <- rep(c("bmi", "whr", "whradjbmi"), each = 3)
types <- rep(c("combined", "females", "males"), 3)
df_names <- paste(traits, "_indexsnps_", sep="")
df_names <- paste(df_names, types, sep="")

file_names <- paste(df_names, ".txt", sep="")

# Read all the tables
dfs <- lapply(file_names, function (x) read.table(paste("Exposures/", x, sep=""), 
                                                  header = T, sep = " ",
                                                  stringsAsFactors = F))
# Name dataframes with their respective titles for easy access
names(dfs) <- df_names

## For the analysis that uses all SNPs, merge the SNPs from combined, females,
# and males analyses to create a single instrument

# To extract dfs by their names, use 'starts with bmi', contains 'whr_', and
# 'whradjbmi' to prevent confounding by overlapping names
traits_grep <- c("^bmi", "^whr_", "whradjbmi")
# Bind dataframes for all BMI SNPs, WHR SNPs, etc.
dfs_all <- lapply(traits_grep, 
                  function (x) do.call("rbind.data.frame", 
                                       dfs[grep(x, names(dfs))]))
# Remove duplicate rows (as would have arisen because SNPs overlap between the
# combined and sex-specific analyses)
dfs_all <- lapply(dfs_all, function(x) x[!duplicated(x), ])

# Name new dfs consistently with previous
traits <- c("bmi", "whr", "whradjbmi")
names(dfs_all) <- paste(traits, "_indexsnps_all", sep="")

# Add the newly collated dataframes to the list of all dataframes
dfs <- c(dfs, dfs_all)

### Format data to build instruments for WHR, BMI, WHRadjBMI ----

# Three ways to weight instruments:
# 1. Female-specific SNPs only with female-specific weights (fem_fem)
# 2. All SNPs with female-specific weights (all_fem)
# 3. All SNPs with combined-sex weights (all_comb)

# Get columns for required for 2-sample MR:
# SNP (rsid), beta, se, effect_allele
# and columns optional for 2-sample MR that are present in the Pulit data
# other_allele, eaf, chr, position, samplesize, and pval

dfs_fem_fem <- lapply(dfs[grep("females", names(dfs))], function (x) 
  x[, c("SNP", "beta.females", "se.females", "A1.combined", "A2.combined",
        "frqA1.females", "Chr.ref.males", "Pos.ref.males", "nmeta.females",
        "pval.females")])
names(dfs_fem_fem) <- paste(traits, "_fem_fem", sep="")

dfs_all_fem <- lapply(dfs[grep("all", names(dfs))], function (x) 
  x[, c("SNP", "beta.females", "se.females", "A1.combined", "A2.combined",
        "frqA1.females", "Chr.ref.males", "Pos.ref.males", "nmeta.females",
        "pval.females")])
names(dfs_all_fem) <- paste(traits, "_all_fem", sep="")

dfs_all_comb <- lapply(dfs[grep("all", names(dfs))], function (x) 
  x[, c("SNP", "beta.combined", "se.combined", "A1.combined", "A2.combined",
        "frqA1.combined", "Chr.ref.males", "Pos.ref.males", "nmeta.combined",
        "pval.combined")])
names(dfs_all_comb) <- paste(traits, "_all_comb", sep="")

dfs <- c(dfs_fem_fem, dfs_all_fem, dfs_all_comb)

# Add column for Phenotype to each df that is useful for 2SMR package
dfs <- mapply(cbind, dfs, "Phenotype" = rep(traits, 3), 
              SIMPLIFY = F)

# Rename columns to standardised names for 2SMR package
changeNames <- function(x) {
  colnames(x) <- c("SNP", "beta", "se", "effect_allele",
                   "other_allele", "eaf", "chr", "position",
                   "samplesize", "pval", "Phenotype")
  return (x)
}

dfs <- lapply(dfs, changeNames)

# Only keep rsid part of SNP, so extract only the portion of SNP before :
getRSID <- function(x) {
  x[, "SNP"] <- sub('\\:.*', '', x[, "SNP"])
  return (x)
}

dfs <- lapply(dfs, getRSID)

# SNPs with negative effect sizes need to be flipped so they have a 
# positive effect size; this involves:
# 1. reversing the sign of beta
# 2. swapping effect allele and other allele
# 3. changing effect allele frequency to 1 - eaf

flipSNPs <- function(x) {
  # get SNPs whose effect sizes are negative
  indexes <- which(x$beta < 0)
  # reverse the sign of beta
  x[indexes, "beta"] <- -x[indexes, "beta"]
  # swap effect allele and other allele
  orig_A1 <- x[indexes, "effect_allele"]
  orig_A2 <- x[indexes, "other_allele"]
  x[indexes, "effect_allele"] <- orig_A2
  x[indexes, "other_allele"] <- orig_A1
  # change eaf to 1 - eaf
  x[indexes, "eaf"] <- 1 - x[indexes, "eaf"]
  
  return (x)
}

dfs <- lapply(dfs, flipSNPs)

# Format data in 2SMR format
# This removes duplicated SNPs
dfs <- lapply(dfs, format_data)

### Clump SNPs ----

# to ensure that the instruments for the exposure are independent:
# European samples from the 1000 genomes project are used to estimate LD 
# between SNPs; amongst those SNPs that have LD R-square above the specified 
# threshold only the SNP with the lowest P-value will be retained
# Automatic cutoffs: window - 10Mb, R2 - 0.001
# Modify to: window - 1Mb, R2 - 0.05 (based on Censin et al. 2019 paper) 

dfs_clumped <- lapply(dfs, function (x) clump_data(x,
                                                   clump_kb = 1000,
                                                   clump_r2 = 0.05))

# List SNPs removed after clumping
removedSNPs <- function (x1, x2) {
  return (x1$SNP[!x1$SNP %in% x2$SNP])
}

# Use sink to divert output from console to text file
sink("Exposures/removed_SNPs.txt")
mapply(removedSNPs, x1 = dfs, x2 = dfs_clumped)
sink()

# ### Read and format data for lean body mass ----
# 
# # Get summary stats from NHGRI GWAS catalog, Zillikens et al. 2017 meta-analysis
# data(gwas_catalog)
# lbm_gwas <- subset(gwas_catalog, grepl("Zillikens", Author))
# 
# # Fill in info for missing SNPs (missing effect_allele and other_allele)
# # using original Zillikens paper
# 
# lbm_gwas[lbm_gwas$SNP == "rs2943656", "effect_allele"] <- "A"
# lbm_gwas[lbm_gwas$SNP == "rs2943656", "other_allele"] <- "G"
# 
# lbm_gwas[lbm_gwas$SNP == "rs4842924", "effect_allele"] <- "T"
# lbm_gwas[lbm_gwas$SNP == "rs4842924", "other_allele"] <- "C"
# 
# lbm_gwas[lbm_gwas$SNP == "rs9936385", "effect_allele"] <- "T"
# lbm_gwas[lbm_gwas$SNP == "rs9936385", "other_allele"] <- "C"
# 
# # Flip SNPs with negative effect sizes
# lbm_gwas <- flipSNPs(lbm_gwas)
# # Get rid of column that says "unit increase" or "unit decrease" as the 
# # SNPs have all been flipped to reflect unit increases
# lbm_gwas <- lbm_gwas[, -23]
# # Format data for 2-Sample MR
# lbm_df <- format_data(lbm_gwas)
# 
# # Add lean body mass data to list of dataframes 
# dfs_clumped$lbm_df <- lbm_df

### Assess instrument strength (F-statistics) ----

calcFStatistic <- function(x) {
  # formula for F-statistic: 1/L * sum over all SNPs(beta.exp^2 / var.exp)
  num <- (x$beta.exposure)^2
  denom <- (x$se.exposure)^2
  fsum <- sum(num/denom)
  fmean <- fsum / length(num)
  return (fmean)
}

# Calculate dilution percentage for each constructed instrument
mean_fstat <- lapply(dfs_clumped, calcFStatistic)
dilution <- lapply(mean_fstat, function (x) 100/x )

# Save F-stats and dilution percentages to file
sink("Exposures/f-stats.txt")
mean_fstat
dilution
sink()

# Method that performs best (highest F-stats) is the sex-specific approach
# that uses female-only SNPs with sex-specific weights
# Continue to use this for all analysis
dfs_winner <- dfs_clumped[grep("fem_fem", names(dfs_clumped))]

# Combine all dataframes into a single "exposures" dataframe for easy analysis
# To do this, need to add "gene.exposure" column to the BMI/WHR/WHRadjBMI dfs
dfs_winner <- mapply(cbind, dfs_winner, "gene.exposure" = NA, 
                     SIMPLIFY = F)

# Retain lean body mass instruments that did not have a weighting strategy
# and a "samplesize.exposure" column to the Lean Body Mass df
dfs_winner$lbm_df <- dfs_clumped$lbm_df
dfs_winner$lbm_df <- within(dfs_winner$lbm_df, {
  samplesize.exposure <- ifelse(exposure == "Lean body mass", 
                                38292, 28330)
})

dfs_winner$lbm_df <- dfs_winner$lbm_df[names(dfs_winner$bmi_fem_fem)]
dfs_winner <- do.call(rbind, dfs_winner)
saveRDS(dfs_winner, "Exposures/dfs_winner.rds")

# Collate outcome data ----

# Extract outcome information for all SNPs present in all instruments; 
# then subset the dataset as needed for each MR

dfs_winner <- readRDS("Exposures/dfs_winner.rds")
# Get list of SNPs to extract information for
all_SNPs <- unique(dfs_winner$SNP)

### Uterine fibroids ----

# Read data from Gallagher et al. 2017 (downloaded from NHGRI GWAS Catalog)
# https://www.ebi.ac.uk/gwas/publications/31649266

uf_outcome <- read_outcome_data(snps = all_SNPs,
                                filename = "Outcomes/Uterine Fibroids/UF_GallagherCS_31649266.stats",
                                sep = "\t", snp_col = "SNP", beta_col = "EFFECT",
                                se_col = "STDERR", 
                                effect_allele_col = "ALLELE1",
                                other_allele_col = "ALLELE0",
                                pval_col = "P", samplesize_col = "N")
# Additional information for phenotype, publication details, etc.
uf_outcome$outcome <- "uterine fibroids"
uf_outcome$data_source.outcome <- "NHGRI GWAS Catalog"

saveRDS(uf_outcome, "Outcomes/Uterine Fibroids/uterine_fibroids.rds")

### Heavy menstrual bleeding ----

# Read data from Gallagher et al. 2017 (downloaded from NHGRI GWAS Catalog)
# https://www.ebi.ac.uk/gwas/publications/31649266

hmb_outcome <- read_outcome_data(snps = all_SNPs,
                                 filename = "Outcomes/Heavy Menstrual Bleeding/GallagherCS_31649266_HeavyBleeding.stats",
                                 sep = "\t", snp_col = "SNP", beta_col = "BETA",
                                 se_col = "SE", 
                                 effect_allele_col = "ALLELE1",
                                 other_allele_col = "ALLELE0",
                                 eaf_col = "A1FREQ",
                                 pval_col = "P_BOLT_LMM")

# Additional information for phenotype, publication details, etc.
hmb_outcome$outcome <- "heavy menstrual bleeding"
hmb_outcome$data_source.outcome <- "NHGRI GWAS Catalog"

saveRDS(hmb_outcome, "Outcomes/Heavy Menstrual Bleeding/hmb.rds")

### UKBB outcomes from SAIGE ----
saige_outcomes <- readRDS("Outcomes/SAIGE/all_outcomes.RDS")

# Keep SAIGE outcomes for: endometriosis, irregular menstrual cycle,,
# female infertility, miscarraiges, and preeclampsia
saige_outcomes <- subset(saige_outcomes, 
                         outcome %in% c("endometriosis", 
                                        "irregular menstrual cycle",
                                        "female infertility", 
                                        "miscarriage and stillbirth",
                                        "preeclampsia"))

# Merge all outcomes datasets
uf_outcome <- readRDS("Outcomes/Uterine Fibroids/uterine_fibroids.rds")
hmb_outcome <- readRDS("Outcomes/Heavy Menstrual Bleeding/hmb.rds")
all_outcomes <- bind_rows(uf_outcome, hmb_outcome)

all_outcomes <- bind_rows(all_outcomes, saige_outcomes)
saveRDS(all_outcomes, "Outcomes/all_outcomes.rds")

# Harmonise data ----

all_outcomes <- readRDS("Outcomes/all_outcomes.rds")
dfs_winner <- readRDS("Exposures/dfs_winner.rds")

harmonised <- harmonise_data(dfs_winner, all_outcomes)

# Rename 'id' columns to later make plots simpler
harmonised$id.exposure <- harmonised$exposure
harmonised$id.outcome <- harmonised$outcome

saveRDS(harmonised, "Results/harmonised.rds")

# Perform MR ----

harmonised <- readRDS("Results/harmonised.rds")

mr_ivw <- mr(harmonised, method_list = "mr_ivw")
# Adjust the obtained p-values for multiple testing
mr_ivw$adj.pval <- p.adjust(mr_ivw$pval, method = "fdr")

# Reorder the dataset to have results ordered by outcome
mr_ivw <- mr_ivw[order(mr_ivw$outcome), ]
sig_results <- subset(mr_ivw, adj.pval <= 0.05)

# For sensitivity, do MR with MR-Egger and weighted median methods
mr_sensitivity <- mr(harmonised, method_list = c("mr_egger_regression",
                                                 "mr_weighted_median"))

all_mr_results <- bind_rows(mr_ivw, mr_sensitivity)
all_mr_results <- all_mr_results[with(all_mr_results, 
                                      order(all_mr_results$outcome,
                                            all_mr_results$exposure)), ]

saveRDS(all_mr_results, "Results/mr_results.rds")
write.table(all_mr_results[, 3:10], "Results/all_results.txt", sep = "\t",
            quote = F, row.names = F)

# Look at significant MR results
unique_combinations <- dim(unique(mr_results[c("outcome", "exposure")]))[1]
p_threshold <- 0.05 / unique_combinations

sig_results <- mr_results[mr_results$pval <= p_threshold, ]

# Sensitivity Analyses ----

# Check for outliers (come back to analyse this if the heterogeneity results
# are off)
uf_snpbysnp <- mr_singlesnp(uf_harmonised)
saveRDS(uf_snpbysnp, "Results/Uterine Fibroids/snp_by_snp.rds")

### Heterogeneity statistics ----

heterogeneity <- mr_heterogeneity(harmonised, 
                                  method_list = c("mr_ivw", 
                                                  "mr_egger_regression"))

# Merge heterogeneity results with MR results
mr_results <- merge(mr_results, heterogeneity[, 3:8], by = c("outcome",
                                                             "exposure",
                                                             "method"),
                    all = T)

# Look for instruments with significant heterogeneity:
sig_heterogeneity <- heterogeneity[heterogeneity$Q_pval <= p_threshold, ]
# There is significant heterogeneity as assessed by the Cochran's Q-statistic
# for some instruments

# As the F-statistics for the instruments were strong and indicated
# little bias (~1.5%), there may be directional pleiotropy
# However, the difference between the Q-statistic for IVW and that for MR-Egger 
# is insignificant for all combinations, suggesting that there is 
# no horizontal pleiotropy

### Horizontal pleiotropy ----

# USE THIS ONCE THERE ARE MULTIPLE OUTCOMES
# uf_hor_pleiotropy <- lapply(uf_harmonised, mr_pleiotropy_test)

hor_pleiotropy <- mr_pleiotropy_test(harmonised)

# There is no significant horizontal pleiotropy - the intercepts for the MR-Egger
# estimate are close to 0 for all combinations
# Merge pleiotropy results with MR results
colnames(hor_pleiotropy)[6:7] <- c("egger_intercept_se", "egger_intercept_pval")
hor_pleiotropy$method <- "MR Egger"
mr_results <- merge(mr_results, hor_pleiotropy[, 3:8], by = c("outcome",
                                                              "exposure",
                                                              "method"),
                    all = T)

## Write all results to table
write.table(mr_results, "Results/all_results.txt", sep = "\t", quote = F, 
            row.names = F)

### Results Plots ----

# Scatter plot to see slopes of various predictors
scatter_plots <- mr_scatter_plot(mr_results, harmonised)
sig_plots <- scatter_plots[paste(sig_results$id.exposure, 
                                 sig_results$id.outcome, sep = ".")]
for (i in 1:length(sig_plots)) {
  filename <- paste("Results/scatter_", 
                    names(sig_plots)[i], sep="")
  ggsave(paste(filename, ".pdf", sep=""), sig_plots[[i]], dpi = 300)
}

# Funnel plot to assess bias / directional pleiotropy
singlesnp <- mr_singlesnp(harmonised)
funnel_plots <- mr_funnel_plot(singlesnp)

# Clearly MR-Egger does not add value to most estimates - intercept is 
# 0 and there is no directional pleiotropy in most of the data;
# this is not true for lean body mass / appendicular lean mass however
# maybe becuase they have so few SNPs to start with
