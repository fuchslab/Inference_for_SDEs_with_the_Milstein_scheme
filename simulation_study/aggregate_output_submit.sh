#!/bin/bash

#SBATCH --output=simulation_study/slurm_messages/output_%x_%j.txt
#SBATCH --error=simulation_study/slurm_messages/error_%x_%j.txt
#SBATCH --job-name=aggregate_output
#SBATCH --partition=icb_cpu
#SBATCH --exclude=ibis-ceph-[002-006,008-019],ibis216-010-[001-004,007,011-012,020-037,051,064,071],ibis216-224-[010-011],icb-rsrv[05-06,08],ibis216-010-[068-070], # remaining:  icb-neu-[001-003]
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=30gb
#SBATCH --nice=0



# #SBATCH --partition=icb_cpu
# #SBATCH --exclude=ibis-ceph-[002-006,008-019],ibis216-010-[001-004,007,011-012,020-037,051,064,071],ibis216-224-[010-011],icb-rsrv[05-06,08] # remaining: ibis216-010-[068-070], icb-neu-[001-003]
# #SBATCH --time=2-00:00:00


# this script is called by simulation_study/aggregate_output_submit_loop.sh
# with argument $1 (dataset/obsFolder) $2 (M) $3 (m) $4 (methodPathUpdate) $5 (approxTransDens) $6 (approxPropDens)


echo dataset: $1
echo M $2
echo m $3
echo methodPathUpdate $4
echo approxTransDens $5
echo approxPropDens $6
echo
echo job name: $SLURM_JOB_NAME
echo
echo job id: $SLURM_JOB_ID
echo array job id: $SLURM_ARRAY_JOB_ID
echo array task id: $SLURM_ARRAY_TASK_ID
echo
echo node name: $SLURMD_NODENAME
echo
echo Current time: $(date)
echo


start=`date +%s`

# start the charliecloud container;
# the user home directory is mounted automatically; however,
# note that the path to the "user home" folder does not contain "/icb" 
    ch-run -b /localscratch:/localscratch/  -b /storage/groups/:/storage/groups  \
    /localscratch/${USER}/r3.6.2_rstan_rmd_rutils/ -- /bin/bash \
     /storage/groups/biostat01/projects/Inference_for_SDEs_with_Milstein/revision_RSOS/simulation_study/aggregate_output_run.sh $1 $2 $3 $4 $5 $6

end=`date +%s`
echo
echo Duration in seconds: $(echo "$end - $start" | bc -l)
echo
echo
echo

if [ ! -d "simulation_study/slurm_messages/${1}/aggregate_output" ] ; then
    mkdir -p simulation_study/slurm_messages/${1}/aggregate_output
fi

# move and rename output and error message files
output_message_file="simulation_study/slurm_messages/${1}/aggregate_output/output_M_${2}_m_${3}_${4}_td_${5}_pd_${6}.txt"
error_message_file="simulation_study/slurm_messages/${1}/aggregate_output/error_M_${2}_m_${3}_${4}_td_${5}_pd_${6}.txt"


cat "simulation_study/slurm_messages/output_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.txt" >> $output_message_file
cat "simulation_study/slurm_messages/error_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.txt"  >> $error_message_file

rm "simulation_study/slurm_messages/output_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.txt"
rm "simulation_study/slurm_messages/error_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.txt"
