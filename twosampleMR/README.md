## Two-sample Mendelian Randomisation to estimate causal effects of obesity traits on female reproductive disorders using publically available GWAS summary statistics

### Files:

Scripts to be run in the following order - 
1. **Scripts in /meta_analysis/** - METAL files to specify inverse-variance weighted fixed effects meta-analysis of UKBIOBANK and FinnGen summary statistics; instructions are specified in *.txt files, and executed from **meta_all_shell.sh**
2. **construct_obesity_instruments.R** - Build BMI, WHR, and WHRadjBMI genetic instruments from index SNPs listed in Pulit et al. (2019), compare different weighting strategies, and assess instrument strength with F-statistics. Use similar script for GIANT-only instruments downloaded from Locke et al. (2015) and Shungin et al. (2015).
3. **get_meta_outcomes.R** - Get meta-analysis outcome results for SNPs in the relevant genetic instruments
4. **obesity_wrh_MR.R** - Perform two-sample MR with constructed instruments and outcomes, including different MR methods and sensitivity tests

For MR with waist circumference, hip circumference, waist-specific WHR, hip-specific WHR, and visceral fat mass instruments, replace step 2. above with: 
1. **construct_fat_instruments.R** - Construct genetic instruments for waist circumference, hip circumference (MRC-IEU database), waist- and hip-specific WHR (Lotta et al.), and predicted visceral fat mass (Karlsson et al.), and calculate instrument strength with F-statistics

For sensitivity analyses with FinnGen-only outcomes, replace step 3. above with: 
1. **get_finngen_R3_outcomes.R** - Get only FinnGen GWAS results for SNPs in the relevant genetic instruments
