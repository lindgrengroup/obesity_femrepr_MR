# Author: Samvida S. Venkatesh
# Date: 27/06/20

library(tidyverse)
theme_set(theme_bw())
library(gridExtra)

## Split by age; read data ----

data <- readRDS("results/summarised_obesity_wrh_phenotypes_210720.rds")
adiposity_bin_start <- as.numeric(str_match(data$adiposity_bin, "[(](.*?)[,]")[, 2])
adiposity_bin_end <- as.numeric(str_match(data$adiposity_bin, "[,](.*?)[]]")[, 2])
data$mean_adiposity <- rowMeans(data.frame(adiposity_bin_start, adiposity_bin_end))

diagnoses <- c("endometriosis", "excessive_menstruation", "infertility", 
               "miscarriage", "PCOS", "preeclampsia.or.eclampsia", 
               "uterine_fibroids")
# Plot

data$lo_ci <- pmax(data$lo_ci, 0)
data$up_ci <- pmin(data$up_ci, 1)

adiposity <- unique(data$adiposity)


# Number of cases and total
data <- data %>% mutate(n_cases = n*mean)

BMI_zoom <- subset(data, data$adiposity == "BMI" & 
                     data$mean_adiposity >= 15 & data$mean_adiposity <= 45)

plots <- list()
for (d in diagnoses) {
  plot_data <- subset(BMI_zoom, BMI_zoom$diagnosis == d)
  plots[[d]] <- ggplot(plot_data, aes(x = mean_adiposity, y = mean)) +
    geom_point(aes(size = n_cases)) +
    geom_line() +
    scale_size(range = c(0.5, 3), guide = F) +
    geom_ribbon(aes(ymin = lo_ci, ymax = up_ci), alpha = 0.2) +
    theme(axis.title = element_blank())
}
tiff("figures/bmi_wrh.tiff", res = 800, units = "cm",
     width = 15, height = 20)
grid.arrange(plots[[1]], plots[[7]], 
             plots[[5]], plots[[2]], 
             plots[[4]], plots[[3]],
             plots[[6]],
             ncol = 2, nrow = 4)
dev.off()

WHR_zoom <- subset(data, data$adiposity == "WHR" & 
                     data$mean_adiposity >= 0.65 & data$mean_adiposity <= 1)
plots <- list()
for (d in diagnoses) {
  plot_data <- subset(WHR_zoom, WHR_zoom$diagnosis == d)
  plots[[d]] <- ggplot(plot_data, aes(x = mean_adiposity, y = mean)) +
    geom_point(aes(size = n_cases)) +
    geom_line() +
    scale_size(range = c(0.5, 3), guide = F) +
    geom_ribbon(aes(ymin = lo_ci, ymax = up_ci), alpha = 0.2) +
    theme(axis.title = element_blank())
}
tiff("figures/whr_wrh.tiff", res = 800, units = "cm",
     width = 15, height = 20)
grid.arrange(plots[[1]], plots[[7]], 
             plots[[5]], plots[[2]], 
             plots[[4]], plots[[3]],
             plots[[6]],
             ncol = 2, nrow = 4)
dev.off()
