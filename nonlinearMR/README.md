## Non-linear Mendelian Randomisation to test for heterogeneity of causal effect of exposure on outcome in different exposure strata

### Files:

Scripts to be run in the following order - 
1. **general_QC.R** - Genotyping and sample quality control in UKBIOBANK data to retain only self-reported and genotype-matched women of white British ancestry that are not related and pass a range of QC checks
2. **calculate_risk_scores.sh** - Construct obesity genetic instruments for each individual that passed QC, with index variants from Pulit et al., using PLINK2 for per-chromosome processing of risk scores. Script here is for BMI, repeat for WHR and WHRadjBMI
3. **combine_risk_scores.R** - Combine PLINK output scores per chromosome into single genetic instrument per individual
4. **collate_nlmr_fields.R** - Gather columns for non-linear MR, including case status (from previous observational phenotype data), obesity status (BMI, WHR, and rank-inverse transformed WHRadjBMI), first 10 principal components, and calculated polygenic risk scores 
5. **perform_nlmr.R** - Calculate variance explained by exposure and perform non-linear MR