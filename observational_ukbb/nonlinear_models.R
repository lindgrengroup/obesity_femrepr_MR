library(mfp)
library(splines)
library(mgcv)
library(tidyverse)
theme_set(theme_bw())

# Read and add data ----

data <- read.table("results/obesity_wrh_phenotypes_passed_qc_210720.txt", 
                   sep = "\t", header = T)

# Make variable for age^2
data$age_sq <- data$age^2

# Ensure that variables for assessment centre and smoking status are factors
data$assessment_centre <- factor(data$assessment_centre)
data$smoking_factor <- factor(ifelse(data$smoking_status == 0, "Never smoked",
                                     ifelse(data$smoking_status == 1 | 
                                              data$smoking_status == 2, "Ever smoked",
                                            NA)))

WHRadjBMImodel <- glm(WHR ~ BMI, data = data, na.action = na.exclude)
data$WHRadjBMI <- residuals(WHRadjBMImodel)
data$WHRadjBMI <- qnorm((rank(data$WHRadjBMI, na.last = "keep") - 0.5) / 
                          sum(!is.na(data$WHRadjBMI)))

diagnoses <- c("endometriosis", "excessive_menstruation", "infertility",
               "miscarriage", "PCOS", "preeclampsia.or.eclampsia", 
               "uterine_fibroids")

# Logistic regression models ----

BMI_logit <- lapply(diagnoses, function (d) {
  m <- bam(formula = paste(d, " ~ BMI + 
        age + age_sq + assessment_centre + smoking_factor", sep = ""),
           data = data, 
           na.action = na.exclude,
           family = binomial(link = "logit"))
  return (m)
})
names(BMI_logit) <- diagnoses

WHR_logit <- lapply(diagnoses, function (d) {
  m <- bam(formula = paste(d, " ~ WHR + 
        age + age_sq + assessment_centre + smoking_factor", sep = ""),
           data = data, 
           na.action = na.exclude,
           family = binomial(link = "logit"))
  return (m)
})
names(WHR_logit) <- diagnoses

WHRadjBMI_logit <- lapply(diagnoses, function (d) {
  m <- bam(formula = paste(d, " ~ WHRadjBMI + 
        age + age_sq + assessment_centre + smoking_factor", sep = ""),
           data = data, 
           na.action = na.exclude,
           family = binomial(link = "logit"))
  return (m)
})
names(WHRadjBMI_logit) <- diagnoses

# Fractional polynomial models ----

BMI_fp <- lapply(diagnoses, function (d) {
  m <- mfp(formula = paste(d, " ~ BMI + 
        age + age_sq + assessment_centre + smoking_factor", sep = ""),
           family = binomial(link = "logit"),
           data = data)
  return (m)
})
names(BMI_fp) <- diagnoses

WHR_fp <- lapply(diagnoses, function (d) {
  m <- mfp(formula = paste(d, " ~ WHR + 
        age + age_sq + assessment_centre + smoking_factor", sep = ""),
           family = binomial(link = "logit"),
           data = data)
  return (m)
})
names(WHR_fp) <- diagnoses

WHRadjBMI_fp <- lapply(diagnoses, function (d) {
  m <- mfp(formula = paste(d, " ~ WHRadjBMI + 
        age + age_sq + assessment_centre + smoking_factor", sep = ""),
           family = binomial(link = "logit"),
           data = data)
  return (m)
})
names(WHRadjBMI_fp) <- diagnoses

fp_models <- lapply(diagnoses, function (d) {
  BMI <- summary(BMI_fp[[d]])
  WHR <- summary(WHR_fp[[d]])
  WHRadjBMI <- summary(WHRadjBMI_fp[[d]])
  return (list(BMI = BMI, WHR = WHR, WHRadjBMI = WHRadjBMI))
})
names(fp_models) <- diagnoses

saveRDS(fp_models, "results/fracpoly_models.rds")

# Generalised additive models ----

BMI_gam <- lapply(diagnoses, function (d) {
  m <- bam(formula = paste(d, " ~ s(BMI) + 
        age + age_sq + assessment_centre + smoking_factor", sep = ""),
           data = data, 
           na.action = na.exclude,
           family = binomial(link = "logit"))
  return (m)
})
names(BMI_gam) <- diagnoses

WHR_gam <- lapply(diagnoses, function (d) {
  m <- bam(formula = paste(d, " ~ s(WHR) + 
        age + age_sq + assessment_centre + smoking_factor", sep = ""),
           data = data, 
           na.action = na.exclude,
           family = binomial(link = "logit"))
  return (m)
})
names(WHR_gam) <- diagnoses

WHRadjBMI_gam <- lapply(diagnoses, function (d) {
  m <- bam(formula = paste(d, " ~ s(WHRadjBMI) + 
        age + age_sq + assessment_centre + smoking_factor", sep = ""),
           data = data, 
           na.action = na.exclude,
           family = binomial(link = "logit"))
  return (m)
})
names(WHRadjBMI_gam) <- diagnoses

# Test GAM model fits vs logistic ----

gam_vs_logit <- lapply(diagnoses, function (d) {
  BMI <- anova(BMI_logit[[d]], BMI_gam[[d]], test = "Chisq")
  WHR <- anova(WHR_logit[[d]], WHR_gam[[d]], test = "Chisq")
  WHRadjBMI <- anova(WHRadjBMI_logit[[d]], WHRadjBMI_gam[[d]], test = "Chisq")
  return (list(BMI = BMI, WHR = WHR, WHRadjBMI = WHRadjBMI))
})

saveRDS(gam_vs_logit, "results/GAM_vs_logit.rds")

# Compute model predictions ----

# As we're interested in the effect of the obesity trait, set the other variables
# to their mean or mode values for factors

# Function to get the mode of a factor vector
getMode <- function (values) {
  # Automatically ignores NA
  t <- data.frame(table(values))
  mode <- t$values[which.max(t$Freq)]
  return (mode)
}

# Create test sets 
testBMIpred <- data.frame(BMI = seq(15, 60, length.out = 1000), 
                          age = mean(data$age, na.rm = T),
                          age_sq = mean(data$age_sq, na.rm = T),
                          assessment_centre = getMode(data$assessment_centre),
                          smoking_factor = getMode(data$smoking_factor))
testWHRpred <- data.frame(WHR = seq(0.6, 1.2, length.out = 1000), 
                          age = mean(data$age, na.rm = T),
                          age_sq = mean(data$age_sq, na.rm = T),
                          assessment_centre = getMode(data$assessment_centre),
                          smoking_factor = getMode(data$smoking_factor))
testWHRadjBMIpred <- data.frame(WHRadjBMI = seq(-4.5, 4.5, length.out = 1000), 
                          age = mean(data$age, na.rm = T),
                          age_sq = mean(data$age_sq, na.rm = T),
                          assessment_centre = getMode(data$assessment_centre),
                          smoking_factor = getMode(data$smoking_factor))

# Predict

BMI_preds <- lapply(diagnoses, function (d) {
  logit_fit <- predict(BMI_logit[[d]], 
                       newdata = testBMIpred, 
                       type = "response", se = T)
  fp_fit <- predict(BMI_fp[[d]], 
                    newdata = testBMIpred, 
                    type = "response", se = T)
  gam_fit <- predict(BMI_gam[[d]], 
                     newdata = testBMIpred, 
                     type = "response", se = T)
  fitted <- data.frame(BMI = testBMIpred$BMI, 
                       fit = logit_fit$fit,
                       lower = logit_fit$fit - 1.96*logit_fit$se.fit,
                       upper = logit_fit$fit + 1.96*logit_fit$se.fit,
                       fit_class = "logistic")
  fitted <- bind_rows(fitted, data.frame(BMI = testBMIpred$BMI, 
                                         fit = fp_fit$fit,
                                         lower = fp_fit$fit - 1.96*fp_fit$se.fit,
                                         upper = fp_fit$fit + 1.96*fp_fit$se.fit,
                                         fit_class = "fracpoly"))
  fitted <- bind_rows(fitted, data.frame(BMI = testBMIpred$BMI, 
                                         fit = gam_fit$fit,
                                         lower = gam_fit$fit - 1.96*gam_fit$se.fit,
                                         upper = gam_fit$fit + 1.96*gam_fit$se.fit,
                                         fit_class = "GAM"))
  return (fitted)
})
names(BMI_preds) <- diagnoses

WHR_preds <- lapply(diagnoses, function (d) {
  logit_fit <- predict(WHR_logit[[d]], 
                       newdata = testWHRpred, 
                       type = "response", se = T)
  fp_fit <- predict(WHR_fp[[d]], 
                    newdata = testWHRpred, 
                    type = "response", se = T)
  gam_fit <- predict(WHR_gam[[d]], 
                     newdata = testWHRpred, 
                     type = "response", se = T)
  fitted <- data.frame(WHR = testWHRpred$WHR, 
                       fit = logit_fit$fit,
                       lower = logit_fit$fit - 1.96*logit_fit$se.fit,
                       upper = logit_fit$fit + 1.96*logit_fit$se.fit,
                       fit_class = "logistic")
  fitted <- bind_rows(fitted, data.frame(WHR = testWHRpred$WHR, 
                                         fit = fp_fit$fit,
                                         lower = fp_fit$fit - 1.96*fp_fit$se.fit,
                                         upper = fp_fit$fit + 1.96*fp_fit$se.fit,
                                         fit_class = "fracpoly"))
  fitted <- bind_rows(fitted, data.frame(WHR = testWHRpred$WHR, 
                                         fit = gam_fit$fit,
                                         lower = gam_fit$fit - 1.96*gam_fit$se.fit,
                                         upper = gam_fit$fit + 1.96*gam_fit$se.fit,
                                         fit_class = "GAM"))
  return (fitted)
})
names(WHR_preds) <- diagnoses

WHRadjBMI_preds <- lapply(diagnoses, function (d) {
  logit_fit <- predict(WHRadjBMI_logit[[d]], 
                       newdata = testWHRadjBMIpred, 
                       type = "response", se = T)
  fp_fit <- predict(WHRadjBMI_fp[[d]], 
                    newdata = testWHRadjBMIpred, 
                    type = "response", se = T)
  gam_fit <- predict(WHRadjBMI_gam[[d]], 
                     newdata = testWHRadjBMIpred, 
                     type = "response", se = T)
  fitted <- data.frame(WHRadjBMI = testWHRadjBMIpred$WHRadjBMI, 
                       fit = logit_fit$fit,
                       lower = logit_fit$fit - 1.96*logit_fit$se.fit,
                       upper = logit_fit$fit + 1.96*logit_fit$se.fit,
                       fit_class = "logistic")
  fitted <- bind_rows(fitted, data.frame(WHRadjBMI = testWHRadjBMIpred$WHRadjBMI, 
                                         fit = fp_fit$fit,
                                         lower = fp_fit$fit - 1.96*fp_fit$se.fit,
                                         upper = fp_fit$fit + 1.96*fp_fit$se.fit,
                                         fit_class = "fracpoly"))
  fitted <- bind_rows(fitted, data.frame(WHRadjBMI = testWHRadjBMIpred$WHRadjBMI, 
                                         fit = gam_fit$fit,
                                         lower = gam_fit$fit - 1.96*gam_fit$se.fit,
                                         upper = gam_fit$fit + 1.96*gam_fit$se.fit,
                                         fit_class = "GAM"))
  return (fitted)
})
names(WHRadjBMI_preds) <- diagnoses

preds <- list(BMI = BMI_preds, WHR = WHR_preds, WHRadjBMI = WHRadjBMI_preds)
saveRDS(preds, "results/nonlinear_model_predictions.rds")
