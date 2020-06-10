#!/bin/bash

# This bash script is executed by 
# `./simulation_study/Stan/aggreagate_stan_output_submit_loop.sh <argument>`
# and one of the arguments 
# [GBM_alpha_1_sigma_2_x0_100 CIR_alpha_1_beta_1_sigma_0.25_x0_3  CIR_alpha_1_beta_1_sigma_2_x0_10]



# number of observations to be used
v_M=(50 20 10)


for M in "${v_M[@]}"
do
  :
    sbatch simulation_study/Stan/aggregate_Stan_output_submit.sh $1 $M 

done
