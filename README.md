# Inference for SDEs with the Milstein scheme
This repository accompanies the article 

Pieschner, S. & Fuchs, C. (2020). Bayesian inference for diffusion processes: using higher-order approximations for transition densities. *Royal Society Open Science*, 7(10), 200270. 

available at [https://royalsocietypublishing.org/doi/10.1098/rsos.200270](https://royalsocietypublishing.org/doi/10.1098/rsos.200270).

**Abstract**

Modelling random dynamical systems in continuous time, diffusion processes are a powerful tool in many areas of science. Model parameters can be estimated from time-discretely observed processes using Markov chain Monte Carlo (MCMC) methods that introduce auxiliary data. These methods typically approximate the transition densities of the process numerically, both for calculating the posterior densities and proposing auxiliary data. Here, the Euler-Maruyama scheme is the standard approximation technique. However, the MCMC method is computationally expensive. Using higher-order approximations may accelerate it, but the specific implementation and benefit remain unclear. Hence, we investigate the utilization and usefulness of higher-order approximations in the example of the Milstein scheme. Our study demonstrates that the MCMC methods based on the Milstein approximation yield good estimation results. However, they are computationally more expensive and can be applied to multidimensional processes only with impractical restrictions. Moreover, the combination of the Milstein approximation and the well-known modified bridge proposal introduces additional numerical challenges.


### *main_functions*
The R-files in the folder *main_functions* implement the parameter estimation methods as described in Pieschner, Fuchs (2020) Section 3.

* `parameter_estimation.R` -- contains the function `estimate_parameters(methodPathUpdate, methodParamUpdate, approxTransDens, approxPropDens, numIterations, m, Y_obs, tau)` that performs the main steps of the estimation algorithm
* `functions_for_parameter_estimation.R` -- contains the general functions used inside the `estimate_parameters()`-function:
    + `update_parameter()`, `update_path()`
    + `rEuler()`, `dEuler()`, `rMilstein()`, `dMilstein()` - implementations of the transition densites based on the two approximation schemes and used for the likelihood density and the left-conditioned proposal densities
    + `rMBEuler()`, `dMBEuler()`  - implementations of the modified bridge proposal densities based on the Euler scheme
    + `rDBMilstein()`, `dDBMilstein()`  - implementations of the diffusion bridge Milstein proposal densities
    + `evaluate_quotient()`, `generateCountingFct()`, `Alg_7_3()`
* `GBM_problem_specific_parameter_and_functions.R`, `CIR_problem_specific_parameter_and_functions.R` -- contain the functions and parameters used inside `estimate_parameters()` that are more specific for our two example, the geometric Brownian motion (GBM) and the Cox-Ingersoll-Ross process (CIR), and our choice of prior and proposal density of the parameter:
    + `len_theta` -- number of parameters (here: 2)
    + `drift_fct()`, `diffusion_fct()`, `diffusion_fct_derivative()`
    + `theta_positive()`-- returns the components of theta with strictly positve range
    + `lambda` -- parameter for the algorithm `Alg_7_3()` to choose the path update interval
    + `rprior_theta()`, `log_dprior_theta()`
    + `propose_theta()` -- random walk proposal density for the parameter vector theta
    + hyperparameters for the prior and proposal density
    + `rMBMilstein()`, `dMBMilstein()` -  functions for the modified bridge proposal densities based on the Milstein scheme are implemented individually as they cannot be easily generalized for all processes

Due to dependencies, the R-files need to be sourced in the following order:

1. `GBM_problem_specific_parameter_and_functions.R` or `CIR_problem_specific_parameter_and_functions.R`
2. `functions_for_parameter_estimation.R`
3. `parameter_estimation.R`

### *simulation_study*
The folder *simulation_study* contains the R-files, input and (aggregated) output files from the simulation study in Section 5 and Appendix D.

* The simulation study was run on a computational grid which is managed with SLURM. The bash scripts `*.sh` were used to submit the individual jobs to the queue of this grid in a loop for the different parameter settings. The computational environment within which the R files were executed is described in `dockerfile_r3.6.2_rstan_rmd_rutils`.
* `observation_generation_GBM.R` (/`observation_generation_CIR.R`) was used to sample the 100 trajectories of the GBM (/CIR) which are saved in one data file `GBM_obs.data` (`GBM_obs.data`) in the folder *GBM_alpha_1_sigma_2_x0_100* (/*CIR_alpha_1_beta_1_sigma_2_x0_10*) along with plots of the trajectories
* `main_simulation_study.R` performs one estimation procedure for the parameters that are passed to it when the script is run and saves a data file of the output to *GBM_alpha_1_sigma_2/output* (/*CIR_alpha_1_beta_1_sigma_0.25/output*)
* `aggregate_output.R`was used to aggregate the results from the individual jobs saved in *GBM_alpha_1_sigma_2/output* (/*CIR_alpha_1_beta_1_sigma_0.25/output*) and save them to the folder *aggregated_output*.
* The folder *Stan* contains all files needed to sample from the true posterior distribution using the Stan software.
* The folder *supplementary_results* contains several .pdf files that give an overview of further results, e.g. results for different numbers M (10, 20, 50) of observed data points.


### *figures\_and\_tables*
The folder *figures\_and\_tables* contains the R-files to generate the figures and tables for the article.
