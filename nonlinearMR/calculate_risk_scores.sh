#!/bin/bash

# This script sets up a task array with one task per operation and uses the step size
# to control how many operations are performed per script run, e.g. to manage the
# turnover time of the tasks. This also makes it a bit easier to re-run a specific
# task than using a step size of one and an unrelated loop counter inside the script

#$ -N BMI_risk_scores
#$ -q short.qc
#$ -P lindgren.prjc
#$ -t 1-22:2
#$ -r y

#$ -cwd -j y
#$ -pe shmem 5

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}
echo SGE_TASK_FIRST=${SGE_TASK_FIRST}, SGE_TASK_LAST=${SGE_TASK_LAST}, SGE_TASK_STEPSIZE=${SGE_TASK_STEPSIZE}

##########################################################################################
#
# Do any one-off set up here

source ~/.bash_profile
mkdir temp_plink_out

##########################################################################################

# Calculate the last task id for this step
this_step_last=$(( SGE_TASK_ID + SGE_TASK_STEPSIZE - 1 ))
if [ "${SGE_TASK_LAST}" -lt "${this_step_last}" ]
then
    this_step_last="${SGE_TASK_LAST}"
fi

# Loop over task ids in this step
while [ "${SGE_TASK_ID}" -le "${this_step_last}" ]
do
    echo `date`: starting work on SGE_TASK_ID=`printenv SGE_TASK_ID`

#   Do per-task processing
	$plink2 --bgen /well/lindgren/UKBIOBANK/DATA/IMPUTATION/ukb_imp_chr${SGE_TASK_ID}_v3.bgen \
	--sample /well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_imp_chr1_v3_s487395.sample \
	--keep results/sample_ids_passed_qc_070720.txt \
	--extract exposures/obesity_instrument_SNPs.txt \
	--score exposures/BMI_instrument_scores.txt 1 2 3 header list-variants ignore-dup-ids cols=denom,scoresums \
	--out temp_plink_out/BMI_GRS_CHR${SGE_TASK_ID}
	rv=$?
# 	Increment SGE_TASK_ID
    export SGE_TASK_ID=$(( SGE_TASK_ID + 1 ))
done

echo `date`: task complete
exit $rv