# Transplacental-antibody-transfer-model
Here we present the first kinetic dynamic model of antibody transfer through the placenta. The components of the model and the results figures are organized as follows:

(0) parameters_Erdogan.m includes the estimated parameter set using CaliPro that was used as the baseline model in the manuscript.  It is taken as input into models in modules 01 and 04.

%%%%%%%%%%%%%%%01-Compartmental Model of Placental Transfer%%%%%%%%%%%%%%%

This folder includes two MATLAB *.m files which call and numerically solve the compartmental IgG transfer model:

(1) diggdt.m includes the set of ordinary differential equations (ODEs), described in Appendix S1 of Erdogan, R. et al

(2) diggdt_driver.m calls and numerically solves the ODEs and plots the output.

%%%%%%%%%%%%%%%02-Parameter Estimation and Sensitivity Analysis%%%%%%%%%%%

This folder includes MATLAB *.m files for modified CaliPro calibration protocol from Joslyn, L. et al (2020) (10.1007/s12195-020-00650-z). Additionally data from Malek, A. et al (1996) (10.1111/j.1600-0897.1996.tb00172.x) is included as CSV files and is called in the calibration protocol.  The files include:

(1) LHS_Call_RE.m - This script is adapted from both Marino, S. et al (2008) (https://doi.org/10.1016/j.jtbi.2008.04.011) and Joslyn, L. et al and performs Latin Hypercube Sampling (LHS) to generate the parameter search space utilized in CaliPro.

(2) RE_Calibration_Driver_pp.m - Driver file for CaliPro. This is the main routine, which takes as input the settings file name (RE_lhs_ode_predator_prey_settings_new_c.m). 

(3) RE_eval_model_performance.m - Function to compare the model runs against the experimental data and determine which of these runs "passed" (according to user specifications) and which runs "failed".

(4) RE_find_ads.m - Function to perform alternative density subtraction to refine parameter space.

(6) RE_lhs_ode_predator_prey_ode_c.m - Contains the ODEs for the compartmental model in a form that accepts parameters sampled using LHS as input.

(7) RE_lhs_ode_predator_prey_settings_new_c.m - This is a settings file that contains initial parameter estimates. The user can edit the number of samples (model runs), timepoints to analyze, initial parameter ranges, and specify the model associated with those parameters.

(8) RE_lhs_ode_run_new_c.m - Function that samples parameter space using a Latin Hypercube Sampling strategy, and runs the model using these parameters.  

(9) PLSR_global_sensitivity_analysis_driver.m - A sample code to run partial least squares discriminant analysis (PLSR) global sensitivity analysis with LHS parameter sets.  This file generates perturbed parameter sets using LHS_Call_RE.m and calls PLSRmain.m, a file which performs PLSR analysis.  The user will need the DolatshahiLab/PLSR-DA_MATLAB repository (DolatshahiLab/PLSR-DA_MATLAB (github.com)) available on GitHub to run this file.

%%%%%%%%%%%%%03-Transwell Assay Simulation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This folder contains an ODE model which simulates IgG transcytosis in HUVEC in an in vitro Transwell assay system.  The files include:

(1) diggdt_transwell.m - Function which contains the ODEs for this model.

(2) parameters_Transwell.m - A file containing parameters used in the Transwell model simulations.

(3) diggdt_transwell_driver.m - Function which loads parameters from parameters_Transwell.m, calls the ODEs in diggdt_transwell.m, numerically solves them, and plots the simulation results.

(4) Erdogan_2023_transwell_data.xlsx - An Excel spreadsheet containing experimental data used to optimize and validate the model, presented in Figures 4 and S6 of the manuscript.

%%%%%%%%%%%%04- Maternal Vaccination and Placental Transfer Combined%%%%

This folder contains two *.m files which analyze the vaccination model.  The files included are:

(1) diggdt_vax.m - Function which contains the ODEs for the vaccine module.  This file is similar to diggdt.m in module 01, but contains equations representing anti-pertussis toxin IgG transfer in addition to bulk maternal IgG transfer.  A new parameter is introduced called p.tvax, which is a user-specified time (given in weeks gestational age) to immunize the mother.

(2) diggdt_vax_driver.m - A file which calls and numerically solves the ODEs in diggdt_vax.m.  The user can specify p.tvax in this script to test different immunization strategies.



