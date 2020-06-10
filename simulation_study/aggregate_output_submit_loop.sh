#!/bin/bash

# This bash script is executed by calling
# `./simulation_study/aggregate_output_submit_loop.sh <argument>`
# and one of the arguments 
# [GBM_alpha_1_sigma_2_x0_100 CIR_alpha_1_beta_1_sigma_0.25_x0_3  CIR_alpha_1_beta_1_sigma_2_x0_10]



# number of observations to be used
v_M=(50 20 10)
# number of imputed points in each inter-observation interval
v_m=(1 2 5)

v_methodPathUpdate=(MB MB MB DBMilstein)
v_approxTransDens=(Euler Milstein Milstein Milstein)
v_approxPropDens=(Euler Euler Milstein Milstein)

v_methodPathUpdate_m1=(leftConditioned leftConditioned)
v_approxTransDens_m1=(Euler Milstein)
v_approxPropDens_m1=(Euler Milstein)


for M in "${v_M[@]}"
do
  :
    for m in "${v_m[@]}"
    do
      :
        if (($m > 1)); then
            for index_method in ${!v_methodPathUpdate[@]};
            do
              :
                methodPathUpdate=${v_methodPathUpdate[index_method]}
                approxTransDens=${v_approxTransDens[index_method]}
                approxPropDens=${v_approxPropDens[index_method]}
                
                sbatch simulation_study/aggregate_output_submit.sh $1 $M $m $methodPathUpdate $approxTransDens $approxPropDens
            done
        else
            for index_method in ${!v_methodPathUpdate_m1[@]};
            do
              :
                methodPathUpdate=${v_methodPathUpdate_m1[index_method]}
                approxTransDens=${v_approxTransDens_m1[index_method]}
                approxPropDens=${v_approxPropDens_m1[index_method]}
                
                sbatch simulation_study/aggregate_output_submit.sh $1 $M $m $methodPathUpdate $approxTransDens $approxPropDens
            done
        fi
    done
done
