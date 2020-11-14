# Author: Samvida S. Venkatesh
# Date: 16/01/2020

PATH = [redacted]

# Load relevant packages
library(TwoSampleMR)   # to perform two-sample MR
library(MRInstruments) # to get GWAS Catalog data for MR

# Waist circumference and hip circumference ----

wc <- TwoSampleMR::extract_instruments(outcomes = "ukb-b-9405",
                                             p1 = 5e-06, clump = T)
wc$exposure <- "Waist circumference"

hc <- TwoSampleMR::extract_instruments(outcomes = "ukb-b-15590",
                                       p1 = 5e-06, clump = T)
hc$exposure <- "Hip circumference"

all_dfs <- list(wc, hc)

### Flip SNPs 

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

all_dfs <- lapply(all_dfs, flipSNPs)

### Assess instrument strength (F-statistics)

calcFStatistic <- function(x) {
  # formula for F-statistic: 1/L * sum over all SNPs(beta.exp^2 / var.exp)
  num <- (x$beta.exposure)^2
  denom <- (x$se.exposure)^2
  fsum <- sum(num/denom, na.rm = T)
  fmean <- fsum / length(num)
  return (fmean)
}

# Calculate dilution percentage for each constructed instrument
lapply(all_dfs, calcFStatistic)

dfs_fat <- bind_rows(all_dfs)

## Waist-sp. and hip-sp. WHR ----

# Construct waist-specific and hip-specific instruments from SNPs in Lotta et al. 2018

waist <- read.table(paste(PATH, "/exposures/fat/waist_specific_whr.txt", sep = ""), 
                    sep = "\t", header = T,
                    stringsAsFactors = F)
waist$Phenotype <- "Waist-specific WHR"
waist$SE <- abs(waist$SE)
hip <- read.table(paste(PATH, "/exposures/fat/hip_specific_whr.txt", sep = ""), 
                  sep = "\t", header = T, 
                  stringsAsFactors = F)
hip$Phenotype <- "Hip-specific WHR"
hip$SE <- abs(hip$SE)

calcFStatistic <- function(x) {
  # formula for F-statistic: 1/L * sum over all SNPs(beta.exp^2 / var.exp)
  num <- (x$Beta)^2
  denom <- (x$SE)^2
  fsum <- sum(num/denom)
  fmean <- fsum / length(num)
  return (fmean)
}

# Calculate dilution percentage for each constructed instrument
mean_fstat <- calcFStatistic(hip)
dilution <- 100/mean_fstat

mean_fstat <- calcFStatistic(waist)
dilution <- 100/mean_fstat

colnames(waist) <- c("SNP", "effect_allele", "other_allele",
                     "eaf", "beta", "se", "pval", "phenotype")
waist <- format_data(waist, type = "exposure")
waist$Phenotype <- "Waist-specific WHR"

colnames(hip) <- c("SNP", "effect_allele", "other_allele",
                   "eaf", "beta", "se", "pval", "phenotype")
hip <- format_data(hip, type = "exposure")
hip$Phenotype <- "Hip-specific WHR"

jama_instruments <- bind_rows(waist, hip)

dfs_fat <- bind_rows(dfs_fat, jama_instruments)

# Visceral fat ----

# Construct visceral fat instruments from SNPs in Karlsson et al.

vfat <- read.table(paste(PATH, "/exposures/fat/vf_karlsson_2019.txt", sep = ""), 
                   sep = "\t", header = T, stringsAsFactors = F)[4:10]

calcFStatistic <- function(x) {
  # formula for F-statistic: 1/L * sum over all SNPs(beta.exp^2 / var.exp)
  num <- (x$BETA)^2
  denom <- (x$SE)^2
  fsum <- sum(num/denom)
  fmean <- fsum / length(num)
  return (fmean)
}

# Calculate dilution percentage for each constructed instrument
mean_fstat <- calcFStatistic(vfat)
dilution <- 100/mean_fstat

colnames(vfat) <- c("SNP", "effect_allele", "eaf",
                    "pval", "beta", "samplesize", "se")
vfat <- format_data(vfat, type = "exposure")
vfat$Phenotype <- "Visceral fat"

dfs_fat <- bind_rows(dfs_fat, vfat)

saveRDS(dfs_fat, paste(PATH, "/exposures/fat/dfs_fat.rds", sep = ""))
