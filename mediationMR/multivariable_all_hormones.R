# Author: Samvida S. Venkatesh
# Date: 21/09/21

library(TwoSampleMR)

lep_FI_instr <- readRDS("exposures/joint_lep_FI_instruments.rds")
lep_ISI_instr <- readRDS("exposures/joint_lep_ISI_instruments.rds")

lep_FI_out <- readRDS("outcomes/lep_FI_outcomes.rds")
lep_ISI_out <- readRDS("outcomes/lep_ISI_outcomes.rds")

exposures <- c("bmi", "whr", "whradjbmi")
outcomes <- c("Endometriosis", "Uterine fibroids", "PCOS",
              "Preeclampsia or eclampsia")

lep_FI_results <- lapply(exposures, function (i) {
  res <- lapply(outcomes, function (o) {
    out <- subset(lep_FI_out[[i]], lep_FI_out[[i]]$outcome == o)
    harmonised <- mv_harmonise_data(lep_FI_instr[[i]], out)
    res <- mv_multiple(harmonised, pval_threshold = 5e-6)#
    res <- as.data.frame(res$result)
    res <- subset(res, res$exposure == i)
    res$exposure <- paste(i, "adj. for leptin and FI")
    return (res)
  })
  res <- bind_rows(res)
  return (res)
})
lep_FI_results <- bind_rows(lep_FI_results)

lep_ISI_results <- lapply(exposures, function (i) {
  res <- lapply(outcomes, function (o) {
    out <- subset(lep_ISI_out[[i]], lep_ISI_out[[i]]$outcome == o)
    harmonised <- mv_harmonise_data(lep_ISI_instr[[i]], out)
    res <- mv_multiple(harmonised, pval_threshold = 5e-6)#
    res <- as.data.frame(res$result)
    res <- subset(res, res$exposure == i)
    res$exposure <- paste(i, "adj. for leptin and ISI")
    return (res)
  })
  res <- bind_rows(res)
  return (res)
})
lep_ISI_results <- bind_rows(lep_ISI_results)

all_res <- bind_rows(lep_FI_results, lep_ISI_results)
all_res$adj_pval <- p.adjust(all_res$pval, method = "fdr")

all_res <- generate_odds_ratios(all_res)
write.table(all_res, "results/joint_mvmr_210921.txt", sep = "\t", quote = F, 
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

calc_for <- read.table("results/tmp.txt", sep = "\t", header = T)

prop_mediated <- deltaMethod(calc_for$EM, calc_for$sigmaEM, 
                             calc_for$MO, calc_for$sigmaMO, 
                             calc_for$EO, calc_for$sigmaEO)
prop_mediated <- bind_cols(calc_for, prop_mediated)

write.table(prop_mediated, "results/two_step_proportion_mediated_030521.txt",
            sep = "\t", quote = F, row.names = F)



