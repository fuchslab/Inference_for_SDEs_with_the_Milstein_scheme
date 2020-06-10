#!/bin/bash

# this script is called by simulation_study/aggregate_output_submit.sh
# with argument $1 (dataset/obsFolder) $2 (M) $3 (m)  $4 (methodPathUpdate) $5 (approxTransDens) $6 (approxPropDens)


echo 
echo output from aggregate_output_run.sh starts here for dataset $1, M $2, m $3, methodPathUpdate $4, approxTransDens $5, approxPropDens $6.
echo

cd /storage/groups/biostat01/projects/Inference_for_SDEs_with_Milstein/revision_RSOS/
  
echo current path: $PWD
echo


Rscript --vanilla simulation_study/aggregate_output.R $1 $2 $3 $4 $5 $6


echo 
echo 
echo output from aggregate_output_run.sh ends here for dataset $1, M $2, m $3, methodPathUpdate $4, approxTransDens $5, approxPropDens $6.