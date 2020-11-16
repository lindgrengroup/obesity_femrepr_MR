## Two-step MR for mediation analysis

### Files:

Scripts to be run in the following order - 
1. **twostep_step1_exposure_mediator.R** - Perform two-sample MR to regress mediators (leptin, fasting insulin, insulin resistance) on exposures (obesity traits) as the first step of two-step mediation MR
2. **twostep_step1_reciprocal.R** - Sensitivity analysis to regress exposures on mediators (reciprocal MR) to confirm direction of causality
3. **mvmr_construct_instruments.R** - Construct genetic instruments for second step of two-step MR (multivariable MR) where genetic instrument for each mediator is adjusted for each obesity trait
4. **mvmr_get_outcomes.R** - Get meta-analysed female reproductive outcome data for SNPs in the MVMR genetic instruments constructed in step 3
5. **twostep_step2_mvmr.R** - Perform multivariable MR (second step of two-step MR) and calculate total, direct, indirect effect sizes and proportion of effect mediated 