#!/bin/sh
echo It is now: $(date)
echo Running on machine: $(hostname)
echo Operating system: $(uname -r)
echo UGE Job ID: $JOB_ID
echo
echo $R_LIBS
cd ~/R/sim_study_CIR/simulation_study/
#cd /mnt/storageGluster/users/susanne.pieschner/R/sim_study/simulation_study

echo $@ $JOB_ID>> grid_output_error_files/combinations_started.txt
echo
echo "The current working directory: $PWD"
echo

echo input parameters:
echo $@

 Rscript main_simulation_study.R $@
echo $@ $JOB_ID>> grid_output_error_files/combinations_finished.txt
echo
echo Now it is: $(date)
echo and Im finished
echo
