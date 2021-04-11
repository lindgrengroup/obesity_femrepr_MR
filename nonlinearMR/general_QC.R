# Author: Samvida S. Venkatesh
# Adapted from: Jenny Censin and Teresa Ferreira
# Date: 07/07/20

# Withdrawn consent ----

remove_withdrawn <- function (data, qc_log_file) {
  
  #Path to UKBB provided list of individuals that have withdrawn consent
  withdrawn <- read.table("/well/lindgren/UKBIOBANK/DATA/QC/w11867_20210201.csv", 
                          header = F)
  
  cleaned <- subset(data, !(data$f.eid %in% withdrawn$V1))
  
  # Report QC metrics
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** Individuals that withdrew consent: ", 
            nrow(data) - nrow(cleaned), "\n",
            "REMAINING, ", nrow(cleaned), "\n", sep = ""))
  sink()
  
  return(cleaned)
  
}

remove_negative_ids <- function (data, qc_log_file) {
  
  cleaned <- subset(data, data$f.eid > 0)
  
  # Report QC metrics
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** Individuals with negative IDs (withdrawn consent): ", 
            nrow(data) - nrow(cleaned), "\n",
            "REMAINING, ", nrow(cleaned), "\n", sep = ""))
  sink()
  
  return(cleaned)
  
}

# Sex ----

qc_sex_mismatch <- function(data, qc_log_file) {
  
  cleaned <- subset(data, 
                    !is.na(data$Submitted.Gender) & !is.na(data$Inferred.Gender) & 
                      data$Submitted.Gender == data$Inferred.Gender)
  
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** Sex mismatch ", 
            nrow(data) - nrow(cleaned), 
            "; REMAINING: ", nrow(cleaned), "\n", sep=""))
  sink()
  
  return(cleaned)	
  
}	

keep_only_females <- function(data, qc_log_file) {
  
  cleaned <- subset(data, data$Inferred.Gender != "M")
  
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** EXCLUDED Males: ", 
            nrow(data) - nrow(cleaned), 
            "; REMAINING: ", nrow(cleaned), "\n", sep=""))
  sink()
  
  return(cleaned)	
  
}

# Ethnicity ----
keep_white_british_ancestry <- function (data, qc_log_file) {
  
  cleaned <- subset(data, data$in.white.British.ancestry.subset == 1)
  
  # Report QC metrics
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** EXCLUDED, Not in white British ancestry subset: ",
            length(which(data$in.white.British.ancestry.subset != 1)), "\n",
            "REMAINING, White British ancestry: ",
            nrow(cleaned), "\n", sep = ""))
  sink()
  
  return (cleaned)
  
}

# Genotyping QC ----

qc_het_miss <- function (data, qc_log_file) {
  
  cleaned <- subset(data, !is.na(data$het.missing.outliers) & 
                      data$het.missing.outliers != 1)	
  
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** EXCLUDED, poor heterozygosity or missingness: ", 
            nrow(data) - nrow(cleaned), 
            "; REMAINING: ", nrow(cleaned), "\n", sep = ""))
  sink()
  
  return(cleaned)
  
}

qc_not_in_phasing  <- function(data, qc_log_file) {
  
  cleaned <- subset(data, !is.na(data$in.Phasing.Input.chr1_22) & 
                      data$in.Phasing.Input.chr1_22 != 0)
  
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** EXCLUDED, not used in autosome phasing: ", 
            nrow(data) - nrow(cleaned), 
            "; REMAINING: ", nrow(cleaned), "\n", sep = ""))
  sink()
  
  return(cleaned)	
  
}		

qc_sex_chr_aneupl  <- function(data, qc_log_file) {
  
  cleaned <- subset(data, !is.na(data$putative.sex.chromosome.aneuploid) & 
                      data$putative.sex.chromosome.aneuploid != 1)
  
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** EXCLUDED, putative sex chr aneuploidy: ", 
            nrow(data) - nrow(cleaned), 
            "; REMAINING: ", nrow(cleaned), "\n", sep = ""))
  sink()
  
  return(cleaned)	
}	

# Relatedness ----

qc_excess_related <- function(data, qc_log_file) {

  cleaned <- subset(data, !is.na(data$excess.relatives) &
                      data$excess.relatives != 1)

  sink(qc_log_file, append = T)
  cat(paste("**FILTER** EXCLUDED, excess relatives (>10 3rd degree relatives): ",
            nrow(data) - nrow(cleaned),
            "; REMAINING: ", nrow(cleaned), "\n", sep = ""))
  sink()

  return(cleaned)

}

qc_related <- function(data, qc_log_file) {

  # Pathway to UKBB list of related individuals
  related <- read.table("/well/lindgren/UKBIOBANK/DATA/QC/ukb1186_rel_s488366.dat",
                        header = T)

  # For each pair of related individuals
  # remove the samples with the highest missingness
  related <- related[related$Kinship > 0.0884 &
                       related$ID1 %in% data$f.eid & related$ID2 %in% data$f.eid, ]

  related$miss1 = data$sample.qc.missing.rate[match(related$ID1, data$f.eid)]
  related$miss2 = data$sample.qc.missing.rate[match(related$ID2, data$f.eid)]
  related$max_miss <- pmax(related$miss1, related$miss2)

  # Remove according to rule above
  related$id_remove <- ifelse(is.na(related$miss1) & is.na(related$miss2),
                              related$ID2,
                              ifelse(is.na(related$miss1), related$ID1,
                                     ifelse(is.na(related$miss2), related$ID2,
                                            ifelse(related$miss1 ==
                                                     related$max_miss, related$ID1,
                                                   ifelse(related$miss2 ==
                                                            related$max_miss,
                                                          related$ID2, "error")))))

  cleaned <- subset(data, !(data$f.eid %in% related$id_remove))

  sink(qc_log_file, append = T)
  cat(paste("**FILTER** Relatedness pairs with errors: ",
            length(which(related$id_remove == "error")), "\n",
            "**FILTER** Individuals excluded because of relatedness: ",
            nrow(data[data$f.eid %in% related$id_remove, ]), "\n",
            "REMAINING NOT RELATED: ", nrow(cleaned), "\n\n", sep = ""))
  sink()

  return(cleaned)

}

qc_kinship_table <- function(data, qc_log_file) {

  cleaned <- subset(data, !is.na(data$excluded.from.kinship.inference) &
                      data$excluded.from.kinship.inference == 0)

  sink(qc_log_file, append = T)
  cat(paste("**FILTER** Excluded from kinship inference: ",
            nrow(data) - nrow(cleaned),
            " ; REMAINING: ", nrow(cleaned), "\n", sep = ""))
  sink()

  return(cleaned)

}

# Other exclusions ----

ukb_recommended_excl <- function (data, qc_log_file) {
  
  # Field: f.22010.0.0, coding: 1 - recommended exclusion
  
  cleaned <- data
  
  remove <- which(cleaned$f.22010.0.0 == 1)
  if (length(remove) > 0) { cleaned <- cleaned[-remove, ] }
  
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** EXCLUDED, recommended UKBIOBANK exclusion: ", 
            nrow(data) - nrow(cleaned), 
            "; REMAINING: ", nrow(cleaned), "\n", sep = ""))
  sink()
  
  return(cleaned)	
  
}

# Prepare data ----

qc_log_file <- "/well/lindgren/UKBIOBANK/samvida/obesity_wrh/nonlinearMR/logs/sample_qc_090421.txt"

# Phenotype file from UKBB
pheno <- read.table("/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv",
                    header = T, sep = ",", na.string = c("NA", "", "."), 
                    stringsAsFactors = F)
colnames(pheno) <- gsub("X", "f.", colnames(pheno))
colnames(pheno)[1] <- "f.eid"

sink(qc_log_file)
cat(paste("SAMPLES IN MAIN PHENOTYPE FILE: ", nrow(pheno), "\n", sep = ""))
sink()

# QC file from UKBB
qc <- read.table("/well/lindgren/UKBIOBANK/DATA/QC/ukb_sqc_v2.txt", header = T, 
                 na.string = c("NA", "", "."), stringsAsFactors = F)
sink(qc_log_file, append = T)
cat(paste("SAMPLES IN UKBB QC FILE: ", nrow(qc), "\n", sep = ""))
cat("\n")
sink()

# fam file corresponding to the QC file provided by the UKBB
fam <- read.table("/well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_cal_chr1_v2_s488363.fam", 
                  header = F)

# Add IDs to QC file
qc$f.eid <- fam[, 1]

# Merge QC file with phenotype file
pheno$Submitted.Gender <- qc$Submitted.Gender[match(pheno$f.eid, qc$f.eid)]
pheno$Inferred.Gender <- qc$Inferred.Gender[match(pheno$f.eid, qc$f.eid)]
pheno$het.missing.outliers <- qc$het.missing.outliers[match(pheno$f.eid, qc$f.eid)]
pheno$excess.relatives <- qc$excess.relatives[match(pheno$f.eid, qc$f.eid)]
pheno$in.Phasing.Input.chr1_22 <- qc$in.Phasing.Input.chr1_22[match(pheno$f.eid, 
                                                                    qc$f.eid)]
pheno$in.white.British.ancestry.subset <- 
  qc$in.white.British.ancestry.subset[match(pheno$f.eid, qc$f.eid)]
pheno$putative.sex.chromosome.aneuploidy <- 
  qc$putative.sex.chromosome.aneuploidy[match(pheno$f.eid, qc$f.eid)]
pheno$sample.qc.missing.rate <- 
  qc$sample.qc.missing.rate[match(pheno$f.eid, qc$f.eid)]
pheno$in.kinship.table <- qc$in.kinship.table[match(pheno$f.eid, qc$f.eid)]
pheno$excluded.from.kinship.inference <- 
  qc$excluded.from.kinship.inference[match(pheno$f.eid, qc$f.eid)]

# Print sample characteristics before QC
sink(qc_log_file, append = T)
cat(paste("TOTAL SAMPLE SIZE: ", nrow(pheno), "\n", sep = ""))
cat(paste("   MEN: ", length(which(pheno$f.31.0.0 == 1)), "\n", sep = ""))
cat(paste("   WOMEN ", length(which(pheno$f.31.0.0 == 0)), "\n", sep = ""))
cat(paste("   GENOTYPED ", length(which(!is.na(pheno$f.22001.0.0))), "\n", sep = ""))
cat(paste("     MEN ", length(which(!is.na(pheno$f.22001.0.0) & 
                                      pheno$f.31.0.0 == 1)), "\n", sep = ""))
cat(paste("     WOMEN ", length(which(!is.na(pheno$f.22001.0.0) & 
                                        pheno$f.31.0.0 == 0)) , "\n", sep = ""))

sink()

# Run QC ----

# Only keep genotyped females
cleaned <- subset(pheno, !is.na(pheno$f.22001.0.0))
cleaned <- keep_only_females(cleaned, qc_log_file)

# Withdrawn
cleaned <- remove_withdrawn(cleaned, qc_log_file)
cleaned <- remove_negative_ids(cleaned, qc_log_file)

# Sex
cleaned <- qc_sex_mismatch(cleaned, qc_log_file)

# Ethnicity
cleaned <- keep_white_british_ancestry(cleaned, qc_log_file)

# Genotyping
cleaned <- qc_het_miss(cleaned, qc_log_file)
cleaned <- qc_not_in_phasing(cleaned, qc_log_file)
cleaned <- qc_sex_chr_aneupl(cleaned, qc_log_file)

# Relatedness
cleaned <- qc_excess_related(cleaned, qc_log_file)
cleaned <- qc_related(cleaned, qc_log_file)
cleaned <- qc_kinship_table(cleaned, qc_log_file)

# Other
cleaned <- ukb_recommended_excl(cleaned, qc_log_file)

# Output ----

# Write file with passed sample IDs
write.table(cleaned$f.eid, 
            "/well/lindgren/UKBIOBANK/samvida/obesity_wrh/nonlinearMR/results/sample_ids_passed_qc_090421.txt", 
            quote = F, row.names = F, sep = "\t")
