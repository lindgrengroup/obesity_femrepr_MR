#!/bin/bash

# Job name
#$ -N FinnGen_UKBB_Meta

# Project name and queue
#$ -P lindgren.prjc 
#$ -q short.qc

# Run the job in the current working directory
#$ -cwd -j y

# Log locations which are relative to the current
# working directory of the submission
###$ -o output.log
###$ -e error.log

# Parallel environment settings (number of nodes / memory)
#$ -pe shmem 3

# Some useful data about the job to help with debugging
echo "------------------------------------------------"
echo "SGE Job ID: $JOB_ID"
echo "SGE Task ID: $SGE_TASK_ID"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"

/apps/well/metal/20110325/metal < ./meta_endometriosis_script.txt
/apps/well/metal/20110325/metal < ./meta_excmens_script.txt
/apps/well/metal/20110325/metal < ./meta_infertility_script.txt
/apps/well/metal/20110325/metal < ./meta_miscarriage_script.txt
/apps/well/metal/20110325/metal < ./meta_pcos_script.txt
/apps/well/metal/20110325/metal < ./meta_pre_or_eclamps_script.txt
/apps/well/metal/20110325/metal < ./meta_uf_script.txt

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0