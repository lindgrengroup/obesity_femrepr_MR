# Author: Samvida S. Venkatesh
# Adapted from: Jenny Censin and Teresa Ferreira
# Date: 17/06/2020

# Withdrawn consent ----

remove_withdrawn <- function (data, qc_log_file) {
  
  #Path to UKBB provided list of individuals that have withdrawn consent
  withdrawn <- read.table("/well/lindgren/UKBIOBANK/DATA/QC/w11867_20200204.csv", 
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

keep_only_females <- function(data, qc_log_file) {
  # Sex is recorded in column 31.0.0, 0 - female, 1 - male
  cleaned <- data[which(data[, "31.0.0"] == 0), ]
  
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** EXCLUDED Males: ", 
            nrow(data) - nrow(cleaned), 
            "; REMAINING: ", nrow(cleaned), "\n", sep=""))
  sink()
  
  return(cleaned)	
  
}

# Ancestry ----
keep_only_white_ancestry <- function (data, qc_log_file) {
  
  # Check for mismatches between assessments and separate European from non-European
  # ancestry
  # Self-reported ethnicity fields - f.21000.0.0, f.21000.1.0, f.21000.2.0
  # Coding for fields:
  # 1 White  1001 - British;  1002 - Irish;  1003 - Any other white background
  # 2 Mixed  2001 - White and Black Caribbean;  2002 - White and Black African;  
  #           2003 - White and Asian;   2004 - Any other mixed background
  # 3 Asian or Asian British 3001 - Indian;  3002 - Pakistani; 3003 - Bangladeshi; 
  #                           3004 - Any other Asian background
  # 4 Black or Black British 4001 - Caribbean; 4002 - African; 
  #                           4003 - Any other Black background
  # 5 Chinese
  # 6 Other ethnic group
  # -1 Do not know
  # -3 Prefer not to answer
  # -10 Exclude
  
  # Get fields
  se_fields <- "f.21000"
  se_fields <- paste(se_fields, c(0, 1, 2), "0", sep = ".")
  
  # Separate individuals of White ancestry from others
  white <- c(1, 1001:1003)
  other <- c(2, 2001:2004, 3, 3001:3004, 4, 4001:4003, 5, 6)
  no_value <- c(0, -1, -3)
  
  data$reported_se <- apply(data[, se_fields], 1, function (x) {
    x[is.na(x)] <- 0
    x <- as.numeric(x)
    ifelse(any(x == -10), "exclude", 
           ifelse(all(is.na(x)), "no_value",
                  ifelse(all(x %in% no_value), "no_value",
                         ifelse(any(x %in% white) & any(x %in% other), "changed_se",
                                ifelse(any(x %in% white), "white", "other")))))
  })
  
  cleaned <- subset(data, data$reported_se == "white")
  
  # Report QC metrics
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** EXCLUDED, -10 in self-reported ethnicity: ",
            length(which(data$reported_se == "exclude")), "\n",
            "**FILTER** EXCLUDED, Has no self-reported ethnicity: ", 
            length(which(data$reported_se == "no_value")), "\n",
            "**FILTER** EXCLUDED, Self-reported mismatch (white v other): ",
            length(which(data$reported_se == "changed_se")), "\n",
            "**FILTER** EXCLUDED, Self-reported other: ", 
            length(which(data$reported_se == "other")), "\n",
            "REMAINING, Self-reported white: ",
            nrow(cleaned), "\n", sep = ""))
  sink()
  
  return (cleaned)

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

qc_log_file <- "/well/lindgren/UKBIOBANK/samvida/obesity_wrh/sample_qc_210720.txt"

# Phenotype file from UKBB
pheno <- read.table("/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv",
                   header = T, sep = ",", na.string = c("NA", "", "."), 
                   stringsAsFactors = F)
colnames(pheno) <- gsub("X", "f.", colnames(pheno))
colnames(pheno)[1] <- "f.eid"

sink(qc_log_file)
cat(paste("SAMPLES IN MAIN PHENOTYPE FILE: ", nrow(pheno), "\n", sep = ""))
sink()

# Print sample characteristics before QC
sink(qc_log_file, append = T)
cat(paste("TOTAL SAMPLE SIZE: ", nrow(pheno), "\n", sep = ""))
cat(paste("   MEN: ", length(which(pheno$f.31.0.0 == 1)), "\n", sep = ""))
cat(paste("   WOMEN ", length(which(pheno$f.31.0.0 == 0)), "\n", sep = ""))
sink()

# Run QC ----

# Withdrawn
cleaned <- remove_withdrawn(pheno, qc_log_file)
cleaned <- remove_negative_ids(cleaned, qc_log_file)

# Sex
cleaned <- keep_only_females(cleaned, qc_log_file)

# Ancestry
cleaned <- keep_only_white_ancestry(cleaned, qc_log_file)

# Other
cleaned <- ukb_recommended_excl(cleaned, qc_log_file)

# Output ----

# Write file with passed sample IDs
write.table(cleaned$f.eid, "/well/lindgren/UKBIOBANK/samvida/obesity_wrh/sample_ids_passed_qc_210720.txt", 
            quote = F, row.names = F, sep = "\t")

# Print sample characteristics after QC
sink(qc_log_file, append = T)
cat(paste("TOTAL SAMPLE SIZE: ", nrow(cleaned), "\n", sep = ""))
cat(paste("   MEN ", length(which(cleaned$f.31.0.0 == 1)), "\n", sep = ""))
cat(paste("   WOMEN ", length(which(cleaned$f.31.0.0 == 0)), "\n", sep = ""))
sink()