## Two-sample Mendelian Randomisation to estimate causal effects of obesity traits on female reproductive disorders using publically available GWAS summary statistics

### Scripts:

Scripts to be run in the following order - 
1. **Scripts in /meta_analysis/** - METAL files to specify inverse-variance weighted fixed effects meta-analysis of UKBIOBANK and FinnGen summary statistics; instructions are specified in *.txt files, and executed from **meta_all_shell.sh**
2. **construct_obesity_instruments.R** - Build BMI, WHR, and WHRadjBMI genetic instruments from index SNPs listed in Pulit et al. (2019), compare different weighting strategies, and assess instrument strength with F-statistics. Use similar script for GIANT-only instruments downloaded from Locke et al. (2015) and Shungin et al. (2015).
3. **get_meta_outcomes.R** - Get meta-analysis outcome results for SNPs in the relevant genetic instruments
4. **obesity_wrh_MR.R** - Perform two-sample MR with constructed instruments and outcomes, including different MR methods and sensitivity tests

For MR with waist circumference, hip circumference, waist-specific WHR, hip-specific WHR, and visceral fat mass instruments, replace step 2. above with: 
1. **construct_fat_instruments.R** - Construct genetic instruments for waist circumference, hip circumference (MRC-IEU database), waist- and hip-specific WHR (Lotta et al.), and predicted visceral fat mass (Karlsson et al.), and calculate instrument strength with F-statistics

For sensitivity analyses with FinnGen-only outcomes, replace step 3. above with: 
1. **get_finngen_R3_outcomes.R** - Get only FinnGen GWAS results for SNPs in the relevant genetic instruments

For sensitivity analyses with multivariable MR for WHR and BMI SNPs in the same model, as compared to WHRadjBMI SNPs - 
1. **Scripts in /whradjbmi_mvmr/** - Scripts to construct instrument with joint exposures, and script to run multivariable MR

### Data:

Harmonised datasets (i.e. alleles aligned across exposures and outcomes) for **main analyses** deposited in /harmonised_datasets/. Main analysis exposures are BMI, WHR, WHRadjBMI (female-specific instruments with female-specific weights, from GIANT-UKB meta-analysis reported by Pulit et al. 2019); other exposures as detailed above. Main analysis outcomes are meta-anlysed across FinnGen, UKB, and/or other largest available European-ancestry GWAS as described in the manuscript S2 Table.  