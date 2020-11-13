# Author: Samvida S. Venkatesh
# Date: 23/06/20

## Read data ----

# Main phenotypes file from UKBB
pheno <- read.table("/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv",
                    header = T, sep = ",", na.string = c("NA", "", "."), 
                    stringsAsFactors = F)
colnames(pheno) <- gsub("X", "f.", colnames(pheno))
colnames(pheno)[1] <- "f.eid"

# IDs to keep that passed QC 
passed_qc <- read.table("/well/lindgren/UKBIOBANK/samvida/obesity_wrh/observational/results/sample_ids_passed_qc_210720.txt",
                        header = T)
colnames(passed_qc) <- "f.eid"

## Extract columns of interest ----

# Only keep samples that passed QC
data <- pheno[pheno$f.eid %in% passed_qc$f.eid, ]

# Columns to extract:
# f.eid - ID
# f.54.0.0 - assessment centre
# f.21003.0.0 - age
# f.3140.0.0 - pregnant (0 - No, 1 - Yes, 2 - Unsure)
# f.20116.0.0 - smoking status (0 - Never, 1 - Previous, 2 - Current, -3 - Prefer not to answer)
# f.21001.0.0 - BMI
# f.48.0.0 - waist circumference
# f.49.0.0 - hip circumference
# f.41202.x.x - main ICD10 diagnosis
# f.41204.x.x - secondary ICD10 diagnosis
# f.40001.x.x - underlying (primary) cause of death ICD10
# f.40002.x.x - contributory (secondary) cause of death ICD10
# f.41203.x.x - main ICD9 diagnosis
# f.41205.x.x - secondary ICD9 diagnosis
# f.20002.x.x - non-cancer illness code, self-reported
# f.2774.x.x - ever had stillbirth, spontaneous miscarriage, or termination (0 - No, 1 - Yes, -1 - Do not know, -3 - Prefer not to answer)

ICD10_cols <- colnames(data)[c(grep("^f.41202.", colnames(data)),
                               grep("^f.41204.", colnames(data)),
                               grep("^f.40001.", colnames(data)),
                               grep("^f.40002.", colnames(data)))]
ICD9_cols <- colnames(data)[c(grep("^f.41203.", colnames(data)),
                              grep("^f.41205.", colnames(data)))]
non_cancer_illness_cols <- colnames(data)[grep("^f.20002.", colnames(data))]
self_reported_miscarr <- colnames(data)[grep("^f.2774.", colnames(data))]

cols_to_extract <- c("f.eid", "f.54.0.0", 
                     "f.21003.0.0", "f.3140.0.0", "f.20116.0.0", 
                     "f.21001.0.0", "f.48.0.0", "f.49.0.0", 
                     ICD10_cols, ICD9_cols, non_cancer_illness_cols, 
                     self_reported_miscarr)
data <- data[, cols_to_extract]
colnames(data) <- c("f.eid", "assessment_centre", 
                    "age", "pregnant", "smoking_status", 
                    "BMI", "wc", "hc", 
                    ICD10_cols, ICD9_cols, non_cancer_illness_cols, 
                    self_reported_miscarr)

## Define phenotypes ----

# WHR
data$WHR <- as.numeric(data$wc) / as.numeric(data$hc)

diagnoses <- c("endometriosis", "excessive_menstruation", "infertility",
               "miscarriage", "PCOS", "preeclampsia or eclampsia", 
               "uterine_fibroids")

# Define ICD10 cases

endometriosis_icd10 <- paste("N80", c("", 0:9), sep = "")
exc_mens_icd10 <- paste("N92", c("", 0:6), sep = "")
infertility_icd10 <- paste("N97", c("", 0, 1, 2, 3, 4, 8, 9), sep = "")
miscarriage_icd10 <- paste("O03", c("", 0:9), sep = "")
PCOS_icd10 <- "E282"
pre_or_eclampsia_icd10 <- c(paste("O14", c("", 0, 1, 2, 9), sep = ""),
                            paste("O15", c("", 0, 1, 2, 9), sep = ""))
uterine_fibroids_icd10 <- paste("D25", c("", 0, 1, 2, 9), sep = "")

icd10_codes <- list(endometriosis_icd10, exc_mens_icd10, infertility_icd10, 
                    miscarriage_icd10, PCOS_icd10, pre_or_eclampsia_icd10,
                    uterine_fibroids_icd10)
names(icd10_codes) <- diagnoses

# Define ICD9 cases

endometriosis_icd9 <- paste("617", c("", 0:9), sep = "")
exc_mens_icd9 <- paste("626", c(2, 3, 4, 5, 6, 8, 9), sep = "")
infertility_icd9 <- paste("628", c("", 0, 1, 2, 3, 4, 8, 9), sep = "")
miscarriage_icd9 <- paste("634", c("", 0:9), sep = "")
PCOS_icd9 <- "2564"
pre_or_eclampsia_icd9 <- paste("642", c(4, 5, 6, 7), sep = "")
uterine_fibroids_icd9 <- paste("218", c("", 9), sep = "")

icd9_codes <- list(endometriosis_icd9, exc_mens_icd9, infertility_icd9, 
                   miscarriage_icd9, PCOS_icd9,
                   pre_or_eclampsia_icd9, uterine_fibroids_icd9)
names(icd9_codes) <- diagnoses

# Define self-reported non-cancer illness cases

endometriosis_nci <- "1402"
exc_mens_nci <- -Inf # Does not exist
infertility_nci <- "1403"
miscarriage_nci <- "1559"
PCOS_nci <- "1350"
pre_or_eclampsia_nci <- "1073"
uterine_fibroids_nci <- c("1351", "1352")

nci_codes <- list(endometriosis_nci, exc_mens_nci, infertility_nci, 
                   miscarriage_nci, PCOS_nci,
                   pre_or_eclampsia_nci, uterine_fibroids_nci)
names(nci_codes) <- diagnoses

# Bring all codes together ----

diagnoses_append <- lapply(diagnoses, function (d) {
  # Look for diagnosis code "d" in relevant columns
  icd10 <- apply(data[, ICD10_cols], 1, 
                 function (x) { any(x %in% icd10_codes[[d]]) })
  
  icd9 <- apply(data[, ICD9_cols], 1, 
                function (x) { any(x %in% icd9_codes[[d]]) })
  
  nci <- apply(data[, non_cancer_illness_cols], 1, 
               function (x) { any(x %in% nci_codes[[d]]) })
  
  # If any of the ICD9, ICD10, or self-reported non-cancer illness codes match, 
  # there is a diagnosis of "d"
  return (icd9 | icd10 | nci)
})
names(diagnoses_append) <- diagnoses

diagnoses_append <- cbind.data.frame(diagnoses_append)
data <- cbind.data.frame(data, diagnoses_append)

# Add cases from primary care data (that were previously calculated)

primarycare_cases <- readRDS("/well/lindgren/UKBIOBANK/samvida/obesity_wrh/observational/results/ukbb_primarycare_wrh_outcomes.rds")

for (d in diagnoses) {
  # Update diagnosis column with primary care results
  pcare <- data$f.eid %in% primarycare_cases[[d]]
  data[, d] <- data[, d] | pcare
}

# Add cases from self reported stillbirth, spontaneous miscarriage, termination

srm <- apply(data[, self_reported_miscarr], 1, function (x) { any(x == 1, na.rm = T)} ) 
data$miscarriage <- data$miscarriage | srm

# Create results data frame ----

res <- data[, c("f.eid", "assessment_centre", 
                "age", "pregnant", "smoking_status", 
                "BMI", "WHR", 
                diagnoses)]
write.table(res, "/well/lindgren/UKBIOBANK/samvida/obesity_wrh/observational/results/obesity_wrh_phenotypes_passed_qc_210720.txt",
            quote = F, row.names = F, sep = "\t")

# Write summary table of number of cases for each phenotype
summary <- colSums(res[, diagnoses])
sink("/well/lindgren/UKBIOBANK/samvida/obesity_wrh/observational/logs/obesity_wrh_phenotypes_passed_qc_210720_summary")
cat(paste("Number of individuals: ", nrow(res), "\n", 
          "Number of cases: ", sep = ""))
cat(paste(names(summary), summary, sep = "- "))
sink()
