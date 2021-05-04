# Author: Samvida S. Venkatesh
# Date: 09/07/20

PATH = [redacted]

library(nlmr)

# Read collated data 
data <- read.table(paste(PATH, "/results/for_nlmr_280720.txt", sep = ""),
                   header = T, sep = "\t", stringsAsFactors = F)

# Create covariates ----

data$age_sq <- data$age^2
# Ensure that variables for assessment centre, genotyping array, and smoking status 
# are factors
data$genotyping.array <- factor(data$genotyping.array)
data$assessment_centre <- factor(data$assessment_centre)
data$smoking_factor <- factor(ifelse(data$smoking_status == 0, "Never smoked",
                                     ifelse(data$smoking_status == 1 | 
                                              data$smoking_status == 2, "Ever smoked",
                                            NA)))

# Remove data with missing values
data <- data[complete.cases(data), ]

covs <- data[, c("age", "age_sq", "genotyping.array", 
                 paste("PC", 1:10, sep = ""),
                 "assessment_centre", "smoking_factor")]

# Trait variance explained by GRS ----

m_bmi_nogrs <- lm(BMI ~ age + age_sq + genotyping.array + 
                   assessment_centre + smoking_factor, data = data)
r1 <- summary(m_bmi_nogrs)$adj.r.squared
m_bmi_grs <- lm(BMI ~ BMI_GRS + age + age_sq + genotyping.array + 
                    assessment_centre + smoking_factor, data = data)
r2 <- summary(m_bmi_grs)$adj.r.squared
# BMI trait variance explained
sink(paste(PATH, "/logs/trait_variance_explained.txt", sep = ""),
     append = T)
cat(paste("R2 BMI model without GRS: ", r1, "\n", sep = ""))
cat(paste("R2 BMI model with GRS: ", r2, "\n", sep = ""))
cat(paste("Variance explained BMI GRS: ", r2 - r1, "\n\n", sep = ""))
sink()

m_whr_nogrs <- lm(WHR ~ age + age_sq + genotyping.array + 
                    assessment_centre + smoking_factor, data = data)
r1 <- summary(m_whr_nogrs)$adj.r.squared
m_whr_grs <- lm(WHR ~ WHR_GRS + age + age_sq + genotyping.array + 
                  assessment_centre + smoking_factor, data = data)
r2 <- summary(m_whr_grs)$adj.r.squared
# WHR trait variance explained
sink(paste(PATH, "/logs/trait_variance_explained.txt", sep = ""),
     append = T)
cat(paste("R2 WHR model without GRS: ", r1, "\n", sep = ""))
cat(paste("R2 WHR model with GRS: ", r2, "\n", sep = ""))
cat(paste("Variance explained WHR GRS: ", r2 - r1, "\n\n", sep = ""))
sink()

m_whradjbmi_nogrs <- lm(WHRadjBMI ~ age + age_sq + genotyping.array + 
                    assessment_centre + smoking_factor, data = data)
r1 <- summary(m_whradjbmi_nogrs)$adj.r.squared
m_whradjbmi_grs <- lm(WHRadjBMI ~ WHRadjBMI_GRS + age + age_sq + genotyping.array + 
                  assessment_centre + smoking_factor, data = data)
r2 <- summary(m_whradjbmi_grs)$adj.r.squared
# WHRadjBMI trait variance explained
sink(pate(PATH, "/logs/trait_variance_explained.txt", sep = ""),
     append = T)
cat(paste("R2 WHRadjBMI model without GRS: ", r1, "\n", sep = ""))
cat(paste("R2 WHRadjBMI model with GRS: ", r2, "\n", sep = ""))
cat(paste("Variance explained WHRadjBMI GRS: ", r2 - r1, "\n\n", sep = ""))
sink()

# Run models ----

diagnoses <- c("endometriosis", "excessive_menstruation", "infertility",
               "miscarriage", "PCOS", "preeclampsia.or.eclampsia", 
               "uterine_fibroids")

fp_bmi <- lapply(diagnoses, function (outcome) {
  y <- as.numeric(data[, outcome])
  x <- data$BMI
  g <- data$BMI_GRS
  
  return (tryCatch(fracpoly_mr(y, x, g, covs, 
                               family = "binomial", q = 100, d = "both", 
                               ci = "model_se", fig = T),
                   error = function (e) NULL))
})

names(fp_bmi) <- diagnoses

saveRDS(fp_bmi, "results/fracpoly.rds")

pl_10 <- lapply(diagnoses, function (outcome) {
  y <- as.numeric(data[, outcome])
  x <- data$BMI
  g <- data$BMI_GRS
  return (tryCatch(piecewise_mr(y, x, g, covs, 
                         family = "binomial", q = 5,
                         fig = T),
                     error = function (e) NULL))
})

saveRDS(pl_10, paste(PATH, "/results/piecewise_q10.rds", sep = ""))