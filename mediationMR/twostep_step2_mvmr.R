# Author: Samvida S. Venkatesh
# Date: 05/08/20

PATH = [redacted]

library(TwoSampleMR)

instruments <- readRDS(paste(PATH, 
                             "/exposures/separate_MVMR_instruments_201016.rds",
                             sep = ""))
outcomes <- readRDS(paste(PATH, "/outcomes/MVMR_outcomes_201016.rds", 
                          sep = ""))

bmi_inst <- grep("^bmi_", names(instruments))
whr_inst <- grep("^whr_", names(instruments))
whradjbmi_inst <- grep("^whradjbmi_", names(instruments))

bmi_results <- lapply(bmi_inst, function (i) {
  relevant_outcomes <- c("uf", "pcos", "pre_or_eclamps")
  res <- lapply(relevant_outcomes, function (o) {
    out <- subset(outcomes[[i]], outcomes[[i]]$outcome == o)
    harmonised <- mv_harmonise_data(instruments[[i]], out)
    res <- mv_multiple(harmonised, pval_threshold = 5e-6)
    return (res)
  })
  names(res) <- relevant_outcomes
  return (res)
})
names(bmi_results) <- names(instruments)[bmi_inst]

whr_results <- lapply(whr_inst, function (i) {
  relevant_outcomes <- c("pcos", "uf")
  res <- lapply(relevant_outcomes, function (o) {
    out <- subset(outcomes[[i]], outcomes[[i]]$outcome == o)
    harmonised <- mv_harmonise_data(instruments[[i]], out)
    res <- mv_multiple(harmonised, pval_threshold = 5e-6)
    return (res)
  })
  return (res)
})

whradjbmi_results <- lapply(whradjbmi_inst, function (i) {
  relevant_outcomes <- c("endometriosis", "uf", "pre_or_eclamps")
  res <- lapply(relevant_outcomes, function (o) {
    out <- subset(outcomes[[i]], outcomes[[i]]$outcome == o)
    harmonised <- mv_harmonise_data(instruments[[i]], out)
    res <- mv_multiple(harmonised, pval_threshold = 5e-6)
    return (res)
  })
  return (res)
})

# Print results to table

all_res <- lapply(list(bmi_results, whr_results, whradjbmi_results), function (j) {
  res <- lapply(j, function (k) {
    res <- lapply(k, function (z) {
      return (z$result)
    })
    res <- bind_rows(res)
  })
  res <- bind_rows(res)
  return (res)
})
all_res <- bind_rows(all_res)

all_res$adj_pval <- p.adjust(all_res$pval, method = "fdr")

all_res <- generate_odds_ratios(all_res)
write.table(all_res, 
            paste(PATH, "/results/mvmr_201016.txt", sep = ""), 
            sep = "\t", quote = F, 
            row.names = F)

# Calculate indirect effects for Two-Step MR ----

# Feed function EM: exposure-mediator estimate, mediator-outcome, and 
# exposure-outcome estimates

deltaMethod <- function(EM, sigmaEM, MO, sigmaMO, EO, sigmaEO) {
  mean_indirect <- EM*MO
  m1 <- eval(D(expression(EM*MO), "EM"))
  m2 <- eval(D(expression(EM*MO), "MO"))
  sigma_indirect <- sqrt((m1^2)*sigmaEM^2 + (m2^2)*sigmaMO^2)
  
  mean_prop_mediated <- mean_indirect / EO
  m3 <- eval(D(expression(mean_indirect / EO), "mean_indirect"))
  m4 <- eval(D(expression(mean_indirect / EO), "EO"))
  sigma_prop_mediated <- sqrt((m3^2)*sigma_indirect^2 + (m4^2)*sigmaEO^2)
  
  res <- list(mean_indirect = mean_indirect, 
              sigma_indirect = sigma_indirect,
              mean_prop_mediated = mean_prop_mediated, 
              sigma_prop_mediated = sigma_prop_mediated,
              lci_prop_mediated = mean_prop_mediated - 1.96*sigma_prop_mediated,
              uci_prop_mediated = mean_prop_mediated + 1.96*sigma_prop_mediated)
  
  return (bind_rows(res))
}

calc_for <- read.table(paste(PATH, "/results/tmp.txt", sep = ""),
                       sep = "\t", header = T)

prop_mediated <- deltaMethod(calc_for$EM, calc_for$sigmaEM, 
                             calc_for$MO, calc_for$sigmaMO, 
                             calc_for$EO, calc_for$sigmaEO)
prop_mediated <- bind_cols(calc_for, prop_mediated)

write.table(prop_mediated, 
            paste(PATH, "/results/two_step_proportion_mediated_201016.txt",
                  sep = ""),
            sep = "\t", quote = F, row.names = F)



