# Author: Samvida S. Venkatesh
# Date: 11/08/20

PATH = [redacted]

library(tidyverse)
theme_set(theme_bw())
library(mrclust)
library(dplyr)
library(pheatmap)

# Read data ----

instruments <- readRDS(paste(PATH, "/data/dfs_winner.rds", sep = ""))
outcomes <- readRDS(paste(PATH, "/data/all_updated_outcomes.rds", sep = ""))

# Harmonise
harmonised <- harmonise_data(instruments, outcomes)

# Run single SNP MR
mr_results <- mr_singlesnp(harmonised)
# Only keep single-SNP results from relevant columns
mr_results <- mr_results[grep("^rs", mr_results$SNP), 
                         c("SNP", "exposure", "outcome",
                           "b", "se")]

# Keep relevant columns for clustering
cluster_data <- harmonised[, c("SNP", "exposure", "outcome",
                             "beta.exposure", "se.exposure",
                             "beta.outcome", "se.outcome")]
cluster_data <- merge(cluster_data, mr_results)

saveRDS(cluster_data, paste(PATH, "/data/SNP_clustering_data.rds", 
                            sep = ""))

# Cluster for each exposure-outcome combination

cluster_data <- readRDS(paste(PATH, "/data/SNP_clustering_data.rds", 
                              sep = ""))

exposures <- c("bmi", "whr", "whradjbmi")
outcomes <- c("Endometriosis", "Heavy Menstrual Bleeding",
              "Infertility", "recurrent miscarriage", 
              "sporadic miscarriage", "PCOS", 
              "Preeclampsia or eclampsia",
              "Uterine fibroids")

cluster_results <- list()
i <- 0

for (exp in exposures) {
  for (out in outcomes) {
    i <- i + 1
    d <- subset(cluster_data, cluster_data$exposure == exp & 
                  cluster_data$outcome == out)
    cluster_results[[i]] <- mr_clust_em(theta = d$b, theta_se = d$se,
                                   bx = d$beta.exposure, 
                                   by = d$beta.outcome,
                                   bxse = d$se.exposure,
                                   byse = d$se.outcome,
                                   obs_names = d$SNP)
  }
}

names(cluster_results) <- paste(rep(exposures, each = 8), 
                                rep(outcomes, times = 3), 
                                sep = " x ")

saveRDS(cluster_results, paste(PATH, "/results/mr_clust_clusters.rds",
                               sep = ""))

# View results - get list of SNPs that belong to a non-junk, non-null cluster
# for each relationship

rel <- names(mr_clust_clusters)
snps <- lapply(1:length(mr_clust_clusters), function (i) {
  m <- mr_clust_clusters[[i]]
  if (!is.null(m)) {
    res <- m$results$best[!m$results$best$cluster_class %in% c("Null", "Junk"), 
                          c("observation", "cluster_class", "probability")]
    keep <- which(res$probability >= 0.8)
    if (length(keep) > 25) {
      res <- res[order(res$probability, decreasing = T), ]
      keep <- 1:25
    }
    res <- res[keep, ]
    res <- res[!is.na(res$observation), ]
    if (dim(res)[1] != 0) {
      res$relationship <- rel[i] 
    }
  }
  return (res)
})

snps <- bind_rows(snps)

saveRDS(snps, paste(PATH, "/results/clustering_top_snps.rds", sep = ""))

# Plot SNP heatmaps ----

# These are the SNP effects in long-form (saved in diseases_clustering.R)
snp_effects <- readRDS(paste(PATH, "/data/causal_effects.rds", sep = ""))
snp_effects$SNP <- as.character(snp_effects$SNP)
# These are the above variants that have been annotated with SNPsnap
annotated_variants <- readRDS(paste(PATH, "/data/annotated_variants.rds",
                                    sep = ""))
annotated_variants$SNP <- as.character(annotated_variants$SNP)
# These are the top causal SNPs generated above
clusters <- readRDS(paste(PATH, "/results/clustering_top_snps.rds", 
                          sep = ""))
clusters$observation <- as.character(clusters$observation)

# For each exposure separately 
exp_list <- c("bmi", "whr", "whradjbmi")

plot_data <- lapply(exp_list, function (exp) {
  
  causal_effects <- subset(snp_effects, snp_effects$exposure == exp) %>% 
    pivot_wider(id_cols = SNP, values_from = b, names_from = outcome)
  # Apply scaling (not mean-centered) to visualise
  d <- causal_effects[, 2:9]
  d <- scale(d, center = F)
  causal_effects[, 2:9] <- d
  
  # Get clusters for given exposure
  grep_string <- paste("^", exp, " x", sep = "")
  cluster_snps <- clusters[grepl(grep_string, clusters$relationship), ]
  cluster_snps$cluster <- sapply(strsplit(cluster_snps$relationship, " x "), 
                                 "[", 2)
  cluster_snps$cluster <- paste(cluster_snps$cluster, 
                                cluster_snps$cluster_class, sep = "_")
  cluster_snps <- cluster_snps[, c("observation", "cluster")]
  colnames(cluster_snps) <- c("SNP", "cluster")
  
  # Merge with gene name
  cluster_snps <- merge(cluster_snps, 
                        distinct(annotated_variants[, c("SNP", "genename")],
                                                    SNP, .keep_all  = T), 
                        by = "SNP", all.x = T)
  cluster_snps$genename <- ifelse(is.na(cluster_snps$genename), 
                                  cluster_snps$SNP, 
                                  cluster_snps$genename)
  
  # Merge with causal effects df per disease
  cluster_snps <- merge(cluster_snps, 
                        causal_effects, by = "SNP", all.x = T)
  
  # Order columns roughly by disease clustering
  cluster_snps <- cluster_snps[, c("SNP", "cluster", "genename",
                                   "Infertility", "Endometriosis",
                                   "Uterine fibroids", "Heavy Menstrual Bleeding", 
                                   "PCOS", "Preeclampsia or eclampsia",
                                   "sporadic miscarriage", "recurrent miscarriage")]
  cluster_snps <- distinct(cluster_snps)
  return (cluster_snps)
})

names(plot_data) <- exp_list

bmi_data <- plot_data$bmi
bmi_data$cluster <- factor(bmi_data$cluster, 
                        levels = c("Endometriosis_1",
                                   "Uterine fibroids_1", 
                                   "Heavy Menstrual Bleeding_1",
                                   "PCOS_1", "Preeclampsia or eclampsia_1"))
# Order by cluster 
bmi_data <- bmi_data[order(bmi_data$cluster), ] 
bmi_annot <- data.frame(cluster = bmi_data$cluster)
annot_gaps <- cumsum(table(bmi_annot$cluster))
plot <- bmi_data[, 4:11]

tiff(paste(PATH, "/figures/mrclust_bmi_clusters.tiff", sep = ""), 
     units = "cm",
     height = 10, width = 15, res = 800)
pheatmap(plot, 
         cluster_rows = F, cluster_cols = F,
         scale = "none", treeheight_col = 0,
         labels_row = bmi_data$genename,
         fontsize = 6,
         show_colnames = F,
         gaps_col = 1:7, gaps_row = annot_gaps)
dev.off()

whr_data <- plot_data$whr
whr_data$cluster <- factor(whr_data$cluster, 
                           levels = c("Endometriosis_1",
                                      "Uterine fibroids_1", 
                                      "Uterine fibroids_2", 
                                      "Heavy Menstrual Bleeding_1",
                                      "Preeclampsia or eclampsia_1"))
# Order by cluster 
whr_data <- whr_data[order(whr_data$cluster), ] 
whr_annot <- data.frame(cluster = whr_data$cluster)
annot_gaps <- cumsum(table(whr_annot$cluster))
plot <- whr_data[, 4:11]

tiff(paste(PATH, "/figures/mrclust_whr_clusters.tiff", sep = ""),
     units = "cm",
     height = 10, width = 15, res = 800)
pheatmap(plot, 
         cluster_rows = F, cluster_cols = F,
         scale = "none", treeheight_col = 0,
         labels_row = whr_data$genename,
         fontsize = 8,
         show_colnames = F,
         gaps_col = 1:7, gaps_row = annot_gaps)
dev.off()


whradjbmi_data <- plot_data$whradjbmi
whradjbmi_data$cluster <- factor(whradjbmi_data$cluster, 
                           levels = c("Endometriosis_1",
                                      "Uterine fibroids_1", 
                                      "Uterine fibroids_2", 
                                      "Heavy Menstrual Bleeding_1",
                                      "Preeclampsia or eclampsia_1"))
# Order by cluster 
whradjbmi_data <- whradjbmi_data[order(whradjbmi_data$cluster), ] 
whradjbmi_annot <- data.frame(cluster = whradjbmi_data$cluster)
annot_gaps <- cumsum(table(whradjbmi_data$cluster))
plot <- whradjbmi_data[, 4:11]

tiff(paste(PATH, "/figures/mrclust_whradjbmi_clusters.tiff", sep = ""),
     units = "cm",
     height = 10, width = 15, res = 800)
pheatmap(plot, 
         cluster_rows = F, cluster_cols = F,
         scale = "none", treeheight_col = 0,
         labels_row = whradjbmi_data$genename,
         fontsize = 8,
         show_colnames = F,
         gaps_col = 1:7, gaps_row = annot_gaps)
dev.off()
