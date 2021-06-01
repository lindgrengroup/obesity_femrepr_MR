# Author: Samvida S. Venkatesh
# Date: 11/08/20

PATH = [redacted]

library(TwoSampleMR)
library(dplyr)
library(umap)
library(tidyverse)
theme_set(theme_bw())

instruments <- readRDS(paste(PATH, "/data/dfs_winner.rds", sep = ""))
outcomes <- readRDS(paste(PATH, "/data/all_updated_outcomes.rds", sep = ""))

# Harmonise
harmonised <- harmonise_data(instruments, outcomes)

# Run single SNP MR
mr_results <- mr_singlesnp(harmonised)
# Only keep single-SNP results
mr_results <- mr_results[grep("^rs", mr_results$SNP), ]
mr_results <- mr_results[, c("SNP", "exposure", "outcome", "b")]

saveRDS(mr_results, paste(PATH, "/data/causal_effects.rds", sep = ""))

# Read data ----

snp_effects <- readRDS(paste(PATH, "/data/causal_effects.rds", sep = ""))

# Look at causal estimate ranges

summ <- snp_effects %>% group_by(exposure, outcome) %>% 
  summarise(min = min(b), median = median(b), 
            mean = mean(b), max = max(b))

# Data is normal within each disease, but standardise the columns
# to see how diseases cluster 
# (otherwise the clusters will be driven by a single disease)

# UMAP for clustering ----

# UMAP for each trait ----

# Replace with BMI, WHR, or WHRadjBMI as appropriate
trait <- "bmi"

d <- subset(snp_effects, snp_effects$exposure == trait)
d <- d[, c("SNP", "outcome", "b")]
d <- d %>% group_by(outcome) %>% 
   pivot_wider(names_from = SNP, values_from = b)
d[is.na(d)] <- 0

d <- data.frame(d, stringsAsFactors = F)
d.labels <- d[, "outcome"]
d <- d[, -1]

# Scale variables to have a variance of 1
d <- t(scale(t(d), center = F))

umap_par <- umap.defaults
umap_par$n_neighbors <- 5
u <- umap(d, config = umap_par)

df <- data.frame(x = u$layout[,1],
                 y = u$layout[,2],
                 Disease = d.labels)

ggplot(df, aes(x, y, label = Disease)) +
   geom_point(size = 3, color = "#c91d59") +
   geom_text(aes(label = Disease), vjust = -0.5, hjust = 1)

tiff("figures/umap_bmi.tiff", width = 7, height = 7, units = "cm", res = 800)
ggplot(df, aes(x, y, label = Disease)) +
   geom_point(color = "#c91d59") +
   geom_hline(yintercept = 0, linetype = "dashed") +
   geom_vline(xintercept = 0, linetype = "dashed") +
   theme(axis.title = element_blank())
dev.off()
