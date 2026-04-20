#!/bin/bash
#
#SBATCH --account=OD-230073
for i in {1..100}
do
JOB_NAME="sim_re${i}"
output="/home/uri003/ReRun_SystematicRecordsinPOdata/out/$JOB_NAME"
sbatch --job-name=$JOB_NAME --output=$output /home/uri003/ReRun_SystematicRecordsinPOdata/job_SurveysinPOdata.slurm $i
done