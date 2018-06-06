#!/bin/sh

# folder containing the file with obvervations
#v_obsFolder=GBM_alpha_1_sigma_2 #(GBM_alpha_1_sigma_2 GBM_alpha_1_sigma_0.25 GBM_alpha_5_sigma_2)
v_obsFolder=(CIR_alpha_1_beta_1_sigma_2 CIR_alpha_1_beta_1_sigma_0.25 CIR_alpha_1_beta_5_sigma_2)
# number of observations to be used
v_M=(50 25)
# number of imputed points in each inter-observation interval
v_m=(1)
# number of iterations
numIterations=1e3

v_methodPathUpdate=(leftConditioned)
v_approxTransDens=(Euler Milstein)
v_approxPropDens=(Euler)


for obsFolder in "${v_obsFolder[@]}"
do
    :
    for index in $(seq 1 10);
    do
        :
        for M in "${v_M[@]}"
        do
            :
            for m in "${v_m[@]}"
            do
                :
                for methodPathUpdate in "${v_methodPathUpdate[@]}"
                do
                    :
                    for approxTransDens in "${v_approxTransDens[@]}"
                    do
                        :
                        for approxPropDens in "${v_approxPropDens[@]}"
                        do
                            :
                            qsub -cwd \
                                -e /home/icb/susanne.pieschner/R/sim_study_CIR/simulation_study/grid_output_error_files/ \
                                -o /home/icb/susanne.pieschner/R/sim_study_CIR/simulation_study/grid_output_error_files/ \
                                -hard -l job_mem=2G -q long_fed25@@susannebenchmark \
                                /home/icb/susanne.pieschner/R/sim_study_CIR/simulation_study/execute_main_file.sh $obsFolder $M $m $numIterations $methodPathUpdate $approxTransDens $approxPropDens $alpha_0 $rho_2 $kappa_0 $nu_0 $gamma_alpha $gamma_sigma $index

                        done
                    done
                done
            done
        done
    done
done
