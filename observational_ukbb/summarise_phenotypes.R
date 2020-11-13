# Author: Samvida S. Venkatesh
# Date: 27/06/20

library(tidyr)
library(dplyr)

PATH = "/well/lindgren/UKBIOBANK/samvida/obesity_wrh/observational/"

## Read and clean data ----
data <- read.table(paste(PATH, "results/obesity_wrh_phenotypes_passed_qc_210720.txt", sep = ""), 
                   sep = "\t", header = T)

# Summarise demographic factors for each diagnosis group ----

# Add a column "no_diagnosis" which is T if no diagnosis
diagnoses <- c("endometriosis", "excessive_menstruation", "infertility",
               "miscarriage", "PCOS", "preeclampsia.or.eclampsia", 
               "uterine_fibroids")
data$no_diagnosis <- rowSums(data[, diagnoses]) == 0

# Gather relevant columns
summ_data <- data[, c("age", "smoking_status", "BMI", "WHR", 
                      "no_diagnosis", diagnoses)]
# Convert smoking status to descriptor
# 0 - Never, 1 - Previous, 2 - Current, -3 - Prefer not to answer
summ_data$smoking_status <- ifelse(summ_data$smoking_status == 0, "Never",
                                   ifelse(summ_data$smoking_status == 1, "Previous",
                                          ifelse(summ_data$smoking_status == 2, 
                                                 "Current", "Not known")))

# Summarise age, BMI, WHR (continuous variables)
grouping_vars <- c(diagnoses, "no_diagnosis")
summ_data <- pivot_longer(summ_data, cols = all_of(grouping_vars),
                          names_to = "diagnosis", values_to = "present")
summ_data <- summ_data[summ_data$present == T, ]

res <- summ_data %>% group_by(diagnosis) %>% 
  summarise(n_cases = n(), 
            mean_age = mean(age, na.rm = T), sd_age = sd(age, na.rm = T),
            mean_BMI = mean(BMI, na.rm = T), sd_BMI = sd(BMI, na.rm = T),
            mean_WHR = mean(WHR, na.rm = T), sd_WHR = sd(WHR, na.rm = T)) %>%
  mutate(LCI_age = mean_age - 1.96*(sd_age/sqrt(n_cases)),
         UCI_age = mean_age + 1.96*(sd_age/sqrt(n_cases)),
         LCI_BMI = mean_BMI - 1.96*(sd_BMI/sqrt(n_cases)),
         UCI_BMI = mean_BMI + 1.96*(sd_BMI/sqrt(n_cases)),
         LCI_WHR = mean_WHR - 1.96*(sd_WHR/sqrt(n_cases)),
         UCI_WHR = mean_WHR + 1.96*(sd_WHR/sqrt(n_cases)))

# Summarise smoking status
res2 <- summ_data %>% group_by(diagnosis, smoking_status) %>% 
  summarise(count = n()) %>% ungroup() %>%
  group_by(diagnosis) %>% 
  mutate(prevalence = count / sum(count))

res2$smoking_status <- paste(res2$smoking_status, "smoker", sep = "_")
res2 <- pivot_wider(res2[, c("diagnosis", "smoking_status",
                             "prevalence")], 
                    names_from = smoking_status, values_from = prevalence)

# Save results
res <- merge(res, res2, by = "diagnosis")
write.table(res, paste(PATH, "logs/summarised_phenotypes_210720.txt", sep = ""),
            sep = "\t", quote = F, row.names = F)

# Prevalence of each diagnosis in BMI and WHR bins ----

# There is one outlier WHR value, so limit WHR to <= 1.5
data <- data[data$WHR <= 1.5, ]
data$BMI_bin <- cut(data$BMI, breaks = 30)
data$WHR_bin <- cut(data$WHR, breaks = 30)

# Pivot longer to summarise for plot
long_data <- pivot_longer(data[, c("BMI_bin", "WHR_bin", diagnoses)],
                          cols = all_of(diagnoses), 
                          names_to = "diagnosis",
                          values_to = "value")
long_data <- subset(long_data, long_data$value %in% c(0, 1))

long_data <- pivot_longer(long_data, 
                          cols = c("BMI_bin", "WHR_bin"),
                          names_to = "adiposity", values_to = "adiposity_bin")
long_data$adiposity <- sub("_bin", "", long_data$adiposity)

# Summarise
plot_adiposity <- long_data %>% group_by(adiposity, adiposity_bin, 
                                             diagnosis) %>% 
  summarise(mean = mean(value), sd = sd(value), n = n()) %>%
  mutate(up_ci = mean + 1.96 * sd/sqrt(n), 
         lo_ci = mean - 1.96 * sd/sqrt(n))

saveRDS(plot_adiposity, paste(PATH, 
                                  "results/summarised_obesity_wrh_phenotypes_210720.rds",
                                  sep = ""))