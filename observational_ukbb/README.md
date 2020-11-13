## Observational associations between obesity traits and female reproductive disorders in UKBIOBANK

### Files:
Scripts to be run in the following order - 
1. **general_QC.R** - UKBIOBANK sample quality control to remove individuals with withdrawn consent, retain only self-identifying females of white British ancestry, and remove UKBB recommended exclusions
2. **icd_to_primarycare.R** - Get primary care read codes for relevant ICD diagnoses codes and extract individuals with female reproductive disease diagnoses in primary care data
3. **get_phenotypes.R** - Collate traits and covariates of interest from UKBIOBANK phenotype file, including age, smoking status, BMI, waist circumference, hip circumference, ICD9 and ICD10 diagnoses, causes of death, self-reported non-cancer illness codes, and self-reported miscarriages
4. **summarise_phenotypes.R** - Summarise demographic and anthropomorphic factors (age, BMI, WHR, etc.) in individuals with and without female reproductive disorders (see Supp. Table 1) 
5. **logit_regression.R** - Perform logistic regression for female reproductive disorders regressed on obesity traits, adjusting for appropriate covariates (see Table 1)
