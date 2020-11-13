# Author: Samvida S. Venkatesh
# Date: 21/07/20

PATH = [redacted]

diagnoses <- c("endometriosis", "excessive_menstruation", "infertility",
               "miscarriage", "PCOS", "preeclampsia or eclampsia", 
               "uterine_fibroids")

# Get ICD 10 codes ----

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


# Get ICD 9 codes ----

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

# Convert ICD codes to read codes for primary care data (v3) ----

v3_icd10_lkps <- read.table(paste(PATH, "/resources/v3_icd10_lkps.txt", sep = ""), sep = "\t", header = T,
                            stringsAsFactors = F)
v3_icd9_lkps <- read.table(paste(PATH, "/resources/v3_icd9_lkps.txt", sep = ""), sep = "\t", header = T,
                           stringsAsFactors = F)

v3_read_codes <- lapply(diagnoses, function (d) {
  res <- v3_icd10_lkps[v3_icd10_lkps$icd10_code %in% icd10_codes[[d]], "read_code"]
  res <- c(res, 
           v3_icd9_lkps[v3_icd9_lkps$icd10_code %in% icd9_codes[[d]], "read_code"])
  return (unique(res))
})
names(v3_read_codes) <- diagnoses

# Convert ICD codes to read codes for primary care data (v2) ----

v2_icd10_lkps <- read.table(paste(PATH, "/resources/v2_icd10_lkps.txt", sep = ""), sep = "\t", header = T,
                            stringsAsFactors = F)
v2_icd9_lkps <- read.table(paste(PATH, "/resources/v2_icd9_lkps.txt", sep = ""), sep = "\t", header = T,
                           stringsAsFactors = F)

v2_read_codes <- lapply(diagnoses, function (d) {
  res <- v2_icd10_lkps[v2_icd10_lkps$icd10_code %in% icd10_codes[[d]], "read_code"]
  res <- c(res, 
           v2_icd9_lkps[v2_icd9_lkps$icd10_code %in% icd9_codes[[d]], "read_code"])
  return (unique(res))
})
names(v2_read_codes) <- diagnoses

saveRDS(v3_read_codes, paste(PATH, "/resources/wrh_primarycare_read3_codes", sep = ""))
saveRDS(v2_read_codes, paste(PATH, "/resources/wrh_primarycare_read2_codes", sep = ""))

# Get f.eid and diagnosis ----

# Read gp_clinical file
gp_clinical <- read.table("/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PRIMARY_CARE/gp_clinical.txt",
                          sep = "\t", header = T, 
                          stringsAsFactors = F)

# Get read codes
v2_read_codes <- readRDS(paste(PATH, "/resources/wrh_primarycare_read3_codes", sep = ""))
v3_read_codes <- readRDS(paste(PATH, "/resources/wrh_primarycare_read3_codes", sep = ""))

diagnoses <- names(v2_read_codes)

# Function to extract whether diagnosis "d" was assigned to a participant
getDiagnosis <- function (d) {
  # v2
  f_eids <- gp_clinical[gp_clinical$read_2 %in% v2_read_codes[[d]], "eid"]
  # v3
  f_eids <- c(f_eids, 
              gp_clinical[gp_clinical$read_3 %in% v3_read_codes[[d]], "eid"])
  # get unique participants with diagnosis of d
  f_eids <- unique(f_eids)
  return (f_eids)
}

# Apply to all diagnoses
wrh_patients <- lapply(diagnoses, function (d) {
  getDiagnosis(d)
})

saveRDS(wrh_patients, paste(PATH, "/results/ukbb_primarycare_wrh_outcomes.rds", sep = ""))
