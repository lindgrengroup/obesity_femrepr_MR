# Author: Samvida S Venkatesh
# Date: 08/06/20

PATH = [redacted]

## Read data ----

# Previously collated phenotypes for observational analysis
pheno <- read.table(paste(PATH, "/results/obesity_wrh_phenotypes_passed_qc_210720.txt", 
                          sep = ""), sep = "\t", header = T, 
                    stringsAsFactors = F)

diagnoses <- c("endometriosis", "excessive_menstruation", "infertility",
               "miscarriage", "PCOS", "preeclampsia.or.eclampsia", 
               "uterine_fibroids")

# QC file from UKBB
qc <- read.table("/well/lindgren/UKBIOBANK/DATA/QC/ukb_sqc_v2.txt", header = T, 
                 na.string = c("NA", "", "."), stringsAsFactors = F)
# fam file corresponding to the QC file provided by the UKBB
fam <- read.table("/well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_cal_chr1_v2_s488363.fam", 
                  header = F)
# Add IDs to QC file
qc$f.eid <- fam[, 1]

# Merge QC file with phenotype file (data on genotyping array and PCs)
pheno$genotyping.array <- qc$genotyping.array[match(pheno$f.eid, qc$f.eid)]
PCs <- paste("PC", 1:40, sep = "")
pheno[, PCs] <- qc[match(pheno$f.eid, qc$f.eid), PCs]

# IDs to keep that passed QC 
passed_qc <- read.table(paste(PATH, "/results/sample_ids_passed_qc_070720.txt", 
                              sep = ""),
                        header = T)
colnames(passed_qc) <- "f.eid"

# Only keep samples that passed QC
data <- pheno[pheno$f.eid %in% passed_qc$f.eid, ]

## Add column for WHRadjBMI ----

# WHRadjBMI (RINT)
WHRadjBMImodel <- lm(WHR ~ BMI, data = data, na.action = na.exclude)
data$WHRadjBMI <- residuals(WHRadjBMImodel)
data$WHRadjBMI <- qnorm((rank(data$WHRadjBMI, na.last = "keep") - 0.5) / 
                          sum(!is.na(data$WHRadjBMI)))

## Append genetic instrument risk scores ----

BMI_GRS <- read.table(paste(PATH, "/results/BMI_GRS.txt", sep = ""),
                      sep = "\t", header = T, 
                      stringsAsFactors = F)
data$BMI_GRS <- BMI_GRS$FINAL_GRS[match(data$f.eid, BMI_GRS$IID)]

WHR_GRS <- read.table(paste(PATH, "/results/WHR_GRS.txt", sep = ""),
                      sep = "\t", header = T, 
                      stringsAsFactors = F)
data$WHR_GRS <- WHR_GRS$FINAL_GRS[match(data$f.eid, WHR_GRS$IID)]

WHRadjBMI_GRS <- read.table(paste(PATH, "/results/WHRadjBMIBMI_GRS.txt", 
                                  sep = ""), sep = "\t", header = T, 
                            stringsAsFactors = F)
data$WHRadjBMI_GRS <- WHRadjBMI_GRS$FINAL_GRS[match(data$f.eid, WHRadjBMI_GRS$IID)]

# Create results data frame ----

res <- data[, c("f.eid", "assessment_centre", "genotyping.array", PCs, 
                "age", "pregnant", "smoking_status", 
                "BMI", "WHR", "WHRadjBMI", 
                "BMI_GRS", "WHR_GRS", "WHRadjBMI_GRS",
                diagnoses)]
write.table(res, paste(PATH, "/results/for_nlmr_280720.txt", sep = ""),
            quote = F, row.names = F, sep = "\t")

# Write summary table of number of cases for each phenotype
summary <- colSums(res[, diagnoses])
sink(paste(PATH, "/logs/for_nlmr_phenotype_count_280720", sep = ""))
cat(paste("Number of individuals: ", nrow(res), "\n", 
          "Number of cases: ", sep = ""))
cat(paste(names(summary), summary, sep = "- "))
sink()