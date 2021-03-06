#!/bin/bash

#SBATCH --output=simulation_study/slurm_messages/output_%A_%a.txt
#SBATCH --error=simulation_study/slurm_messages/error_%A_%a.txt
#SBATCH --job-name=Stan_Milstein_revision
#SBATCH --partition=icb_cpu
#SBATCH --exclude=ibis-ceph-[002-006,008-019],ibis216-010-[001-004,007,011-012,020-037,051,064],ibis216-224-[010-011],icb-neu-[001-003],icb-rsrv[05-06,08] # remaining: ibis216-010-[068-070,071]
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=50gb
#SBATCH --nice=0
#SBATCH --array=1,2,4,16,17,20,24,26,33,39,41,44,45,46,49,50,54,57,62,68,71,73,76,78,81,82,83,84,85,87,89,91,92,93,94,96,97,98,99,100


# This bash script can be submitted to SLURM with 
# `sbatch simulation_study/Stan/stan_submit.sh <argument>`
# and one of the arguments 
# [GBM_alpha_1_sigma_2_x0_100 CIR_alpha_1_beta_1_sigma_0.25_x0_3  CIR_alpha_1_beta_1_sigma_2_x0_10]

echo job name: $SLURM_JOB_NAME
echo dataset: $1
echo
echo job id: $SLURM_JOB_ID
echo array job id: $SLURM_ARRAY_JOB_ID
echo array task id: $SLURM_ARRAY_TASK_ID
echo
echo node name: $SLURMD_NODENAME
echo
echo Current time: $(date)
echo
echo home directory: $HOME

if [ ! -d "/localscratch/${USER}/r3.6.2_rstan_rmd_rutils" ] ; then
    # create a local folder for your own user
    # -p: don't report error if it already exists
    mkdir -p /localscratch/${USER}/tmp
    
    # extract the container image from a tarball
    ch-tar2dir ~/charliecloud_containers/r3.6.2_rstan_rmd_rutils.tar.gz /localscratch/${USER}
fi


    start=`date +%s`

    # start the charliecloud container;
    # the user home directory is mounted automatically; however,
    # note that the path to the "user home" folder does not contain "/icb" 
    ch-run -b /localscratch:/localscratch/  -b /storage/groups/:/storage/groups  \
    /localscratch/${USER}/r3.6.2_rstan_rmd_rutils/ -- /bin/bash \
     /storage/groups/biostat01/projects/Inference_for_SDEs_with_Milstein/revision_RSOS/simulation_study/Stan/stan_run.sh $1

    end=`date +%s`
    echo
    echo Duration in seconds: $(echo "$end - $start" | bc -l)
    echo
    echo
    echo
    

if [ ! -d "simulation_study/slurm_messages/Stan/${1}" ] ; then
    mkdir -p simulation_study/slurm_messages/Stan/${1}
fi

# move and rename output and error message files
output_message_file="simulation_study/slurm_messages/Stan/${1}/output_$SLURM_ARRAY_TASK_ID.txt"
error_message_file="simulation_study/slurm_messages/Stan/${1}/error_$SLURM_ARRAY_TASK_ID.txt"


cat "simulation_study/slurm_messages/output_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt" >> $output_message_file
cat "simulation_study/slurm_messages/error_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt"  >> $error_message_file

rm "simulation_study/slurm_messages/output_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt"
rm "simulation_study/slurm_messages/error_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt"



