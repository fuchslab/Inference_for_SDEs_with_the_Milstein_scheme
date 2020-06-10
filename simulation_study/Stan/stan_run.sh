#!/bin/bash

# this script is called by ~/R/Inference_for_SDEs_with_Milstein/simulation_study/Stan/stan_submit.sh
# with argument $1 (dataset/obsFolder) 

# redefine path for a tmp directory (the default tmp directory is mounted as non-executable on the new cluster)
export TMP=/localscratch/${USER}/tmp && export TMPDIR=$TMP && export TEMP=$TMP

echo 
echo output from run.sh starts here for dataset $1 and trajectory $SLURM_ARRAY_TASK_ID.
echo
echo TMP: $TMP

cd /storage/groups/biostat01/projects/Inference_for_SDEs_with_Milstein/revision_RSOS/

if [ ! -d "simulation_study/${1}/output/Stan/stanfit_objects" ] ; then
    mkdir -p simulation_study/${1}/output/Stan/stanfit_objects
    mkdir -p simulation_study/${1}/output/Stan/optimizing_results
fi


echo current path: $PWD
echo

# number of iterations
numIterations=5e5

# number of observations to be used
v_M=(50) #(10 20 50)


for M in "${v_M[@]}"
do
  :           
    echo Current time: $(date)   
            
    output_Rdata=simulation_study/${1}/output/Stan/stanfit_objects/stanfit_object_M_${M}_path_$SLURM_ARRAY_TASK_ID.rds
    # check whether output already exists
    if ls $output_Rdata 1> /dev/null 2>&1 ; then
        echo $(ls $output_Rdata)
        echo already exists
        echo
    else
        echo
        echo Run Stan sampling for
        echo $1 $M $numIterations $SLURM_ARRAY_TASK_ID
        echo 
        
        Rscript --vanilla simulation_study/Stan/stan_sampling.R $1 $M $numIterations $SLURM_ARRAY_TASK_ID 
        
        echo
    fi
done



echo 
echo output from stan_run.sh ends here for obsFolder $1 and trajectory $SLURM_ARRAY_TASK_ID.
echo
