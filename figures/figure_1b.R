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
