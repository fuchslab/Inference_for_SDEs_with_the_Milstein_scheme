#!/bin/bash

# this script is called by simulation_study/Stan/aggregate_Stan_output_submit.sh
# with argument $1 (dataset/obsFolder) $2 (M) 

echo 
echo output from aggregate_output_run.sh starts here for dataset $1, M $2.
echo

cd /storage/groups/biostat01/projects/Inference_for_SDEs_with_Milstein/revision_RSOS/
  
echo current path: $PWD
echo


Rscript --vanilla simulation_study/Stan/aggregate_Stan_output.R $1 $2 


echo 
echo 
echo output from aggregate_output_run.sh ends here for dataset $1, M $2.