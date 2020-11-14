# Author: Samvida S. Venkatesh
# Date: 16/01/2020

PATH = [redacted]

# Load relevant packages
library(tidyverse)     # for data manipulation
library(TwoSampleMR)   # to perform two-sample MR
library(MRInstruments) # to get GWAS Catalog data for MR

# Construct Instruments ----
### Read exposure data for WHR, BMI, WHRadjBMI ----

# Read Pulit et al. 2018 data
# Name all the data by trait and whether the SNPs are from combined, male, or 
# female analyses

traits <- rep(c("bmi", "whr", "whradjbmi"), each = 3)
types <- rep(c("combined", "females", "males"), 3)
df_names <- paste(traits, "_indexsnps_", sep="")
df_names <- paste(df_names, types, sep="")

file_names <- paste(df_names, ".txt", sep="")

# Read all the tables
dfs <- lapply(file_names, function (x) read.table(paste(PATH, x, sep = ""), 
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
sink(paste(PATH, "/logs/obesity_instruments_removed_SNPs_240720.txt", 
           sep = ""))
mapply(removedSNPs, x1 = dfs, x2 = dfs_clumped)
sink()

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
sink(paste(PATH, "/obesity_fstats.txt", sep = ""))
mean_fstat
dilution
sink()

# Method that performs best (highest F-stats) is the sex-specific approach
# that uses female-only SNPs with sex-specific weights
# Continue to use this for all analysis
dfs_winner <- dfs_clumped[grep("fem_fem", names(dfs_clumped))]

# Sensitivity analysis: all SNPs with combined scores
dfs_sensitivity <- dfs_clumped[grep("all_comb", names(dfs_clumped))]

# Combine all dataframes into a single "exposures" dataframe for easy analysis
# To do this, need to add "gene.exposure" column to the BMI/WHR/WHRadjBMI dfs
dfs_winner <- mapply(cbind, dfs_winner, "gene.exposure" = NA, 
                     SIMPLIFY = F)

dfs_sensitivity <- mapply(cbind, dfs_clumped, "gene.exposure" = NA, SIMPLIFY = F)

dfs_winner <- do.call(rbind, dfs_winner)
saveRDS(dfs_winner, paste(PATH, "dfs_winner.rds", sep = ""))


dfs_sensitivity <- do.call(rbind, dfs_sensitivity)
saveRDS(dfs_sensitivity, paste(PATH, "dfs_sensitivity.rds", sep = ""))