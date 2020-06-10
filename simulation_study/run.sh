#!/bin/bash

# this script is called by simulation_study/submit.sh
# with argument $1 (dataset/obsFolder) 

# redefine path for a tmp directory (the default tmp directory is mounted as non-executable on the new cluster)
export TMP=/localscratch/${USER}/tmp && export TMPDIR=$TMP && export TEMP=$TMP

echo 
echo output from run.sh starts here for dataset $1 and trajectory $SLURM_ARRAY_TASK_ID.
echo
echo TMP: $TMP

cd /storage/groups/biostat01/projects/Inference_for_SDEs_with_Milstein/revision_RSOS/

echo current path: $PWD
echo


# computational time before timeout
compTime=3600

# number of iterations
numIterations=1e5

# number of observations to be used
v_M=(50 20 10) # 50
# number of imputed points in each inter-observation interval
v_m=(2 5)

v_methodPathUpdate=(MB MB MB DBMilstein)
v_approxTransDens=(Euler Milstein Milstein Milstein)
v_approxPropDens=(Euler Euler Milstein Milstein)

v_methodPathUpdate_m1=(leftConditioned leftConditioned)
v_approxTransDens_m1=(Euler Milstein)
v_approxPropDens_m1=(Euler Milstein)



for M in "${v_M[@]}"
do
  :
    # m > 1
    for m in "${v_m[@]}"
    do
      :
        for index_method in $(seq 0 3);
        do
          :
            methodPathUpdate=${v_methodPathUpdate[index_method]}
            approxTransDens=${v_approxTransDens[index_method]}
            approxPropDens=${v_approxPropDens[index_method]}
        
            echo Current time: $(date)   
            
            output_Rdata=simulation_study/${1}/output/${1:0:2}*_M_${M}_m_${m}_path_${methodPathUpdate}_td_${approxTransDens}_pd_${approxPropDens}_$SLURM_ARRAY_TASK_ID.data
            # check whether output already exists
            if ls $output_Rdata 1> /dev/null 2>&1 ; then
                echo $(ls $output_Rdata)
                echo already exists
                echo
            else
                echo
                echo Run algorithm for
                echo $1 $M $m $numIterations $methodPathUpdate $approxTransDens $approxPropDens  $SLURM_ARRAY_TASK_ID $compTime
                echo 
                
                timeout $(expr $compTime + 1000) Rscript --vanilla simulation_study/main_simulation_study.R $1 $M $m $numIterations $methodPathUpdate $approxTransDens $approxPropDens  $SLURM_ARRAY_TASK_ID $compTime
                
                echo
            fi
        done
    done
    
    # m = 1
    m=1
    for index_method in $(seq 0 1);
    do
        :
        methodPathUpdate=${v_methodPathUpdate_m1[index_method]}
        approxTransDens=${v_approxTransDens_m1[index_method]}
        approxPropDens=${v_approxPropDens_m1[index_method]}
        
       echo Current time: $(date)
        
        output_Rdata=simulation_study/${1}/output/${1:0:2}*_M_${M}_m_${m}_path_${methodPathUpdate}_td_${approxTransDens}_pd_${approxPropDens}_$SLURM_ARRAY_TASK_ID.data
        # check whether output already exists
        if ls $output_Rdata 1> /dev/null 2>&1 ; then
            echo $(ls $output_Rdata)
            echo already exists
            echo
        else
            echo
            echo Run algorithm for
            echo $1 $M $m $numIterations $methodPathUpdate $approxTransDens $approxPropDens  $SLURM_ARRAY_TASK_ID $compTime
            echo 
            
            timeout $(expr $compTime + 1000) Rscript --vanilla simulation_study/main_simulation_study.R $1 $M $m $numIterations $methodPathUpdate $approxTransDens $approxPropDens  $SLURM_ARRAY_TASK_ID $compTime
            
            echo
        fi
    done
done

echo 
echo output from run.sh ends here for obsFolder ${1} and trajectory $SLURM_ARRAY_TASK_ID.
echo
