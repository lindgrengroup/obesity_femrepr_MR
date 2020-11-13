# Author: Samvida S. Venkatesh
# Date: 30/06/20

library(dplyr)

PATH = [redacted]

## Read and clean data ----
data <- read.table(paste(PATH, "/results/obesity_wrh_phenotypes_passed_qc_210720.txt", sep = ""), 
                   sep = "\t", header = T)

# Adjust traits for covariates ----

# Remove extreme outliers: 
# BMI <= 15 and BMI >= 60
data <- data[data$BMI > 15 & data$BMI < 60, ]
# WHR <= 0.6 and WHR >= 1.2
data <- data[data$WHR > 0.6 & data$WHR < 1.2, ]

# Make variable for age^2
data$age_sq <- data$age^2

# Ensure that variables for assessment centre and smoking status are factors
data$assessment_centre <- factor(data$assessment_centre)
data$smoking_factor <- factor(ifelse(data$smoking_status == 0, "Never smoked",
                                     ifelse(data$smoking_status == 1 | 
                                              data$smoking_status == 2, "Ever smoked",
                                            NA)))

# BMI adjusted:
BMImodel <- glm(BMI ~ age + age_sq + assessment_centre + smoking_factor,
               data = data,
               na.action = na.exclude)
data$adjBMI <- residuals(BMImodel)
data$adjBMI <- qnorm((rank(data$adjBMI, na.last = "keep") - 0.5) / 
                       sum(!is.na(data$adjBMI)))

# WHR adjusted:
WHRmodel <- glm(WHR ~ age + age_sq + assessment_centre + smoking_factor,
               data = data,
               na.action = na.exclude)
data$adjWHR <- residuals(WHRmodel)
data$adjWHR <- qnorm((rank(data$adjWHR, na.last = "keep") - 0.5) / 
                       sum(!is.na(data$adjWHR)))

# WHRadjBMI:
WHRadjBMImodel <- glm(WHR ~ BMI + age + age_sq + assessment_centre + smoking_factor,
               data = data,
               na.action = na.exclude)
data$WHRadjBMI <- residuals(WHRadjBMImodel)
data$WHRadjBMI <- qnorm((rank(data$WHRadjBMI, na.last = "keep") - 0.5) / 
                       sum(!is.na(data$WHRadjBMI)))

# Run logistic regression ----

diagnoses <- c("endometriosis", "excessive_menstruation", "infertility",
               "miscarriage", "PCOS", "preeclampsia.or.eclampsia", 
               "uterine_fibroids")

BMI_model_fits <- lapply(diagnoses, function (d) {
  m <- glm(formula = paste(d, " ~ adjBMI", sep = ""),
      data = data, 
      na.action = na.exclude,
      family = binomial(link = "logit"))
  beta <- summary(m)$coefficients[, 1][2]
  sd <- summary(m)$coefficients[, 2][2]
  pval <- summary(m)$coefficients[, 4][2]
  d <- data.frame(obesity_trait = "BMI", diagnosis = d, 
                  OR = exp(beta), 
                  LCI = exp(beta) - 1.96*sd, UCI = exp(beta) + 1.96*sd,
                  pval = pval)
  return (d)
})
BMI_model_fits <- bind_rows(BMI_model_fits)

WHR_model_fits <- lapply(diagnoses, function (d) {
  m <- glm(formula = paste(d, " ~ adjWHR", sep = ""),
      data = data, 
      na.action = na.exclude,
      family = binomial(link = "logit"))
  beta <- summary(m)$coefficients[, 1][2]
  sd <- summary(m)$coefficients[, 2][2]
  pval <- summary(m)$coefficients[, 4][2]
  d <- data.frame(obesity_trait = "WHR", diagnosis = d, 
                  OR = exp(beta), 
                  LCI = exp(beta) - 1.96*sd, UCI = exp(beta) + 1.96*sd,
                  pval = pval)
  return (d)
})
WHR_model_fits <- bind_rows(WHR_model_fits)

WHRadjBMI_model_fits <- lapply(diagnoses, function (d) {
  m <- glm(formula = paste(d, " ~ WHRadjBMI", sep = ""),
      data = data, 
      na.action = na.exclude,
      family = binomial(link = "logit"))
  beta <- summary(m)$coefficients[, 1][2]
  sd <- summary(m)$coefficients[, 2][2]
  pval <- summary(m)$coefficients[, 4][2]
  d <- data.frame(obesity_trait = "WHRadjBMI", diagnosis = d, 
                  OR = exp(beta), 
                  LCI = exp(beta) - 1.96*sd, UCI = exp(beta) + 1.96*sd,
                  pval = pval)
  return (d)
})
WHRadjBMI_model_fits <- bind_rows(WHRadjBMI_model_fits)

all <- bind_rows(BMI_model_fits, WHR_model_fits, WHRadjBMI_model_fits)

all$adj_pval <- p.adjust(all$pval, method = "fdr")

write.table(all, paste(PATH, "/results/logistic_regression_results_210720.txt", sep = ""),
            sep = "\t", quote = F, row.names = F)
