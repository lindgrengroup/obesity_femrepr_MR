BMI_quantiles <- readRDS("results/obesity_trait_quantiles.rds")[[1]]

tiff("figures/obesity_quantiles.tiff", units = "cm",
     height = 8, width = 8, res = 800)
ggplot(BMI_quantiles, aes(x = value, col = q, fill = q)) +
  geom_density(alpha = 0.2, position = "identity") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_reverse(lim = c(50, 15)) +
  coord_flip() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) 
dev.off()

piecewise_q5 <- readRDS("results/piecewise_q5.rds")

diagnoses <- c("endometriosis", "excessive_menstruation", "infertility",
               "miscarriage", "PCOS", "preeclampsia.or.eclampsia", 
               "uterine_fibroids")
obesity_traits <- c("BMI", "WHR", "WHRadjBMI")

names(piecewise_q5) <- diagnoses

piecewise_q5 <- lapply(piecewise_q5, function (x) {
  names(x) <- obesity_traits
  return (x)
})

d <- data.frame(piecewise_q5[["uterine_fibroids"]][["BMI"]]$lace)
d$q <- factor(1:5)
d$OR <- exp(d$beta)
d$LCI_OR <- exp(d$lci)
d$UCI_OR <- exp(d$uci)
d$significant <- as.factor(ifelse(d$pval < 0.05, 1, 2))

tiff("figures/BMI_UF_lace.tiff", units = "cm",
     height = 7.6, width = 8.2, res = 800)
ggplot(d, aes(x = q, y = OR, ymin = LCI_OR, ymax = UCI_OR)) +
  geom_point(aes(col = q, alpha = significant),
             size = 3) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_errorbar(aes(ymin = LCI_OR, ymax = UCI_OR, 
                    col = q, linetype = significant,
                    alpha = significant), cex = 1) +
  scale_y_log10() +
  scale_x_discrete(lim = rev(levels(d$q))) +
  scale_color_brewer(palette = "Dark2") +
  scale_alpha_discrete(range = c(1, 0.5)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.title = element_blank(),
        legend.position = "none") +
  coord_flip()
dev.off()

for (d in diagnoses) { 
  for (o in obesity_traits) {
    p <- piecewise_q5[[d]][[o]]$figure
    p$theme <- theme_bw() +
      theme(axis.title = element_blank(),
            axis.text = element_blank())
    
    tiff(paste("figures/", d, "_", o, ".tiff", sep = ""), units = "cm",
         height = 12, width = 15, res = 800)
    print (p)
    dev.off()
  }
}
