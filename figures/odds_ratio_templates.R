# Author: Samvida S. Venkatesh
# Date: 25/06/20

library(tidyverse)
theme_set(theme_bw())
library(paletteer)
library(gridExtra)

# Figure 1 - method robustness ----

pub <- read.table("Results/obesity_meta_wrh_mr_results_200930.txt", header = T,
                  sep = "\t")
waisthip <- read.table("Results/waisthip_wrh_results_200930.txt", header = T,
                       sep = "\t")
vfat <- read.table("Results/vfat_wrh_results_200930.txt", header = T,
                   sep = "\t")
vfat$exposure <- "Visceral fat"

data <- bind_rows(pub, vfat, waisthip)

plot_data <- data[, c("outcome", "exposure", 
                      colnames(data)[grep("^OR", colnames(data))], 
                      colnames(data)[grep("^LCI", colnames(data))],
                      colnames(data)[grep("^UCI", colnames(data))],
                      colnames(data)[grep("^pval", colnames(data))])]

# Format data
plot_data <- pivot_longer(plot_data, cols = -c("outcome", "exposure"), 
                          names_to = "Value_Method", values_to = "number")
plot_data <- separate(plot_data, Value_Method, into = c("Value", "Method"), 
                      sep = "_")
plot_data <- pivot_wider(plot_data, id_cols = -c("Value", "number"),
                         names_from = Value, values_from = number)
plot_data$Method <- ifelse(plot_data$Method == "weighted", "Weighted Median", 
                           plot_data$Method)
plot_data$exposure <- ifelse(plot_data$exposure == "bmi", "BMI",
                             ifelse(plot_data$exposure == "whr", "WHR",
                                    ifelse(plot_data$exposure == "whradjbmi", "WHRadjBMI",
                                           plot_data$exposure)))
plot_data$significant <- ifelse(plot_data$pval < 0.05, 1, 2)
plot_data$exposure <- factor(plot_data$exposure, levels = c("BMI", 
                                                            "Waist circumference",
                                                            "Hip circumference",
                                                            "WHR",
                                                            "Waist-specific WHR",
                                                            "Hip-specific WHR",
                                                            "WHRadjBMI", 
                                                            "Visceral fat"))
plot_data$significant <- factor(plot_data$significant)
plot_data$Method <- factor(plot_data$Method, levels = c("Weighted Median", "Egger",
                                                        "IVW"))
p <- list()
for (o in unique(plot_data$outcome)) {
  
  tmp <- subset(plot_data, plot_data$outcome == o)
  
  p[[o]] <- ggplot(tmp, aes(x = Method, y = OR, ymin = LCI, ymax = UCI)) +
    facet_wrap(~exposure, strip.position = "left", nrow = 3, scales = "free") +
    geom_point(aes(col = Method, alpha = significant, shape = Method),
               size = 2) +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_errorbar(aes(ymin = LCI, ymax = UCI, col = Method, linetype = significant,
                      alpha = significant), cex = 0.7) +
    scale_y_log10() +
    scale_color_paletteer_d("nord::frost") +
    scale_alpha_discrete(range = c(1, 0.5)) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(face = "bold"),
          axis.title = element_blank(),
          strip.text.y = element_blank(),
          legend.position = "none") +
    coord_flip()
}

# For recurrent miscarriage, as there are no significant associations, need to 
# specify alpha and linetype

o <- "recurrent miscarriage"
tmp <- subset(plot_data, plot_data$outcome == o)

p[[o]] <- ggplot(tmp, aes(x = Method, y = OR, ymin = LCI, ymax = UCI)) +
  facet_wrap(~exposure, strip.position = "left", nrow = 3, scales = "free") +
  geom_point(aes(col = Method, shape = Method), alpha = 0.5, 
             size = 2) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_errorbar(aes(ymin = LCI, ymax = UCI, col = Method), 
                alpha = 0.5, linetype = 2, cex = 0.7) +
  scale_y_log10() +
  scale_color_paletteer_d("nord::frost") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.title = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none") +
  coord_flip()

tiff("Figures/robust_methods/all_methods_set1.tiff", units = "cm", 
     height = 14, width = 15, res = 800)
grid.arrange(p[[1]], p[[2]], p[[4]], p[[8]], ncol = 2)
dev.off()

# Figure 2 - instrument robustness ---- 

female <- read.table("Results/obesity_meta_wrh_mr_results_200930.txt", header = T,
                     sep = "\t")
female$instrument <- "female"
combined <- read.table("Results/combined_obesity_wrh_results_200930.txt", header = T,
                       sep = "\t")
combined$instrument <- "combined sexes"

dat <- bind_rows(female, combined)

plot_data <- dat[, c("outcome", "exposure", "instrument")]
plot_data$OR <- ifelse(dat$method == "IVW", dat$OR_IVW, 
                       ifelse(dat$method == "Egger", dat$OR_Egger,
                              dat$OR_weighted_median))
plot_data$LCI <- ifelse(dat$method == "IVW", dat$LCI_IVW, 
                        ifelse(dat$method == "Egger", dat$LCI_Egger,
                               dat$LCI_weighted_median))
plot_data$UCI <- ifelse(dat$method == "IVW", dat$UCI_IVW, 
                        ifelse(dat$method == "Egger", dat$UCI_Egger,
                               dat$UCI_weighted_median))
plot_data$pval <- ifelse(dat$method == "IVW", dat$pval_IVW, 
                         ifelse(dat$method == "Egger", dat$pval_Egger,
                                dat$pval_weighted_median))
plot_data$exposure <- ifelse(plot_data$exposure == "bmi", "BMI",
                             ifelse(plot_data$exposure == "whr", "WHR",
                                    "WHRadjBMI"))
plot_data$exposure <- factor(plot_data$exposure, levels = c("WHRadjBMI", "WHR",
                                                            "BMI"))
plot_data$outcome <- factor(plot_data$outcome, 
                            levels = c("Endometriosis", "Uterine fibroids",
                                       "PCOS", "Heavy Menstrual Bleeding",
                                       "sporadic miscarriage", "recurrent miscarriage",
                                       "Infertility", "Preeclampsia or eclampsia"))
plot_data$significant <- ifelse(plot_data$pval < 0.05, 1, 2)

plot_data$significant <- factor(plot_data$significant)
plot_data$yaxis <- paste(plot_data$exposure, plot_data$instrument)
plot_data$yaxis <- factor(plot_data$yaxis, 
                          levels = c("WHRadjBMI combined sexes", "WHRadjBMI female",
                                     "WHR combined sexes", "WHR female", 
                                     "BMI combined sexes", "BMI female"))

tiff("Figures/robust_instruments/all_outcomes.tiff", units = "cm", 
     height = 10, width = 15, res = 800)
ggplot(plot_data, aes(x = yaxis, y = OR, ymin = LCI, ymax = UCI)) +
  facet_wrap(~outcome, strip.position = "left", nrow = 3, scales = "free") +
  geom_point(aes(col = exposure, alpha = significant, shape = instrument),
             size = 2) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_errorbar(aes(ymin = LCI, ymax = UCI, col = exposure, linetype = significant,
                    alpha = significant), cex = 0.7) +
  scale_y_log10() +
  scale_color_paletteer_d("ggsci::default_aaas") +
  scale_alpha_discrete(range = c(1, 0.5)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.title = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none") +
  coord_flip()
dev.off()

# Figure 3 - outcome robustness ----

meta <- read.table("Results/obesity_meta_wrh_mr_results_200930.txt", header = T,
                     sep = "\t", stringsAsFactors = F)
meta$outcome_type <- "meta"
meta <- subset(meta, !meta$outcome %in% c("sporadic miscarriage", 
                                          "recurrent miscarriage"))
finngen <- read.table("Results/obesity_finngen_wrh_mr_results_200728.txt", header = T,
                       sep = "\t", stringsAsFactors = F)
finngen$outcome_type <- "finngen"
finngen <- subset(finngen, finngen$outcome != "Miscarriage")
finngen$outcome <- ifelse(finngen$outcome == "Excessive or frequent menstruation",
                          "Heavy Menstrual Bleeding", finngen$outcome)

dat <- bind_rows(meta, finngen)

plot_data <- dat[, c("outcome", "exposure", "outcome_type")]
plot_data$OR <- ifelse(dat$method == "IVW", dat$OR_IVW, 
                       ifelse(dat$method == "Egger", dat$OR_Egger,
                              dat$OR_weighted_median))
plot_data$LCI <- ifelse(dat$method == "IVW", dat$LCI_IVW, 
                        ifelse(dat$method == "Egger", dat$LCI_Egger,
                               dat$LCI_weighted_median))
plot_data$UCI <- ifelse(dat$method == "IVW", dat$UCI_IVW, 
                        ifelse(dat$method == "Egger", dat$UCI_Egger,
                               dat$UCI_weighted_median))
plot_data$pval <- ifelse(dat$method == "IVW", dat$pval_IVW, 
                         ifelse(dat$method == "Egger", dat$pval_Egger,
                                dat$pval_weighted_median))
plot_data$exposure <- ifelse(plot_data$exposure == "bmi", "BMI",
                             ifelse(plot_data$exposure == "whr", "WHR",
                                    "WHRadjBMI"))
plot_data$exposure <- factor(plot_data$exposure, levels = c("WHRadjBMI",
                                                            "WHR",
                                                            "BMI"))
plot_data$outcome <- factor(plot_data$outcome, 
                            levels = c("Endometriosis", "Uterine fibroids",
                                       "PCOS", "Heavy Menstrual Bleeding",
                                       "Infertility", "Preeclampsia or eclampsia"))
plot_data$significant <- ifelse(plot_data$pval < 0.05, 1, 2)

plot_data$significant <- factor(plot_data$significant)
plot_data$yaxis <- paste(plot_data$exposure, plot_data$outcome_type)
plot_data$yaxis <- factor(plot_data$yaxis, 
                          levels = c("WHRadjBMI finngen", "WHRadjBMI meta",
                                     "WHR finngen", "WHR meta", 
                                     "BMI finngen", "BMI meta"))

tiff("Figures/robust_outcomes/all_outcomes.tiff", units = "cm", 
     height = 7, width = 13, res = 800)
ggplot(plot_data, aes(x = yaxis, y = OR, ymin = LCI, ymax = UCI)) +
  facet_wrap(~outcome, strip.position = "left", nrow = 2, scales = "free") +
  geom_point(aes(col = exposure, alpha = significant, shape = outcome_type),
             size = 2) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_errorbar(aes(ymin = LCI, ymax = UCI, col = exposure, linetype = significant,
                    alpha = significant), cex = 0.7) +
  scale_y_log10() +
  scale_color_paletteer_d("ggsci::default_aaas") +
  scale_alpha_discrete(range = c(1, 0.5)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.title = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none") +
  coord_flip()
dev.off()

