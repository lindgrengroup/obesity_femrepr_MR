# Author: Samvida S. Venkatesh
# Date: 21/09/21

library(TwoSampleMR)

instruments <- readRDS("exposures/mvmr_whr_and_bmi_instr.rds")
outcomes <- readRDS("Outcomes/outcomes_200921.rds")

outcomes <- subset(outcomes, outcomes$SNP %in% instruments$SNP)

# Harmonise and perform MVMR ----

harmonised <- mv_harmonise_data(instruments, outcomes)
res <- mv_multiple(harmonised, pval_threshold = 5e-6)
res <- res$result
whradjbmi_res <- subset(res, res$exposure == "whr")

whradjbmi_res$adj_pval <- p.adjust(whradjbmi_res$pval, 
                                   method = "fdr")

whradjbmi_res <- generate_odds_ratios(all_res)
write.table(whradjbmi_res, 
            "Results/mvmr_whradjbmi_results_210921.txt",
            sep = "\t", quote = F, 
            row.names = F)
