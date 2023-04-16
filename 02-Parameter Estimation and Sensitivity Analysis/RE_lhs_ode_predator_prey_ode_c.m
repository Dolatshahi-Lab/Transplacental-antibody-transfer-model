% This script sourced from:  http://malthus.micro.med.umich.edu/CaliPro/
% Joslyn, L. et al. CaliPro:  A Calibration Protocol That Utilizes Parameter
% Density Estimation to Explore Parameter Space and Calibrate Complex 
% Biological Models.  Cellular and Molecular Bioengineering (2021).
% doi.org/10.1007/s12195-020-00650-z

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% And adapted by:  Erdogan, R. A quantitative mechanistic model reveals
% key determinants of maternal-fetal IgG transfer with implications for 
% prenatal immunization (2023).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The ODEs for a simple predator/prey model with 2 equations and 4 parameters.
% This model is from page 6, section 3.1 of "A Methodology For Performing
% Global Uncertainty And Sensitivity Analysis In Systems Biology", Marino, et.
% al., Journal of Theoretical Biology, 2008-09-07, doi:
% 10.1016/j.jtbi.2008.04.011.

% From the above cited paper, with y(1) = Q and y(2) = P.

% It has two state variables (Q, P) and several parameters. Q represents the
% density of prey, P represents the density of predators, alpha is the intrinsic
% rate of prey population increase, beta is the predation rate coefficient, sigma is
% the predator mortality rate, and delta is the reproduction rate of predators per
% prey consumed.

% This file can be used as a template for creating a new ODE model. To do so
% copy this file and edit it, replacing the predator/prey parameters and
% equations with those for the new model.
%
% When creating a new model, you will also need to copy and edit file
% lhs_ode_predator_prey_settings_new.m to create a settings file for the new
% model, to define the model parameters values, initial conditions, etc.

function dy = RE_lhs_ode_predator_prey_ode_c(t, y, params)

neq = size(y,1);
dy = zeros(neq,1);

% ODE equations
% dy(1) = -params.k_up*y(1)/params.v_m; %maternal blood
% dy(2) = -params.k_up*y(2)/params.v_m; %maternal blood
% dy(3) = -params.k_up*y(3)/params.v_m; %maternal blood
% dy(4) = -params.k_up*y(4)/params.v_m; %maternal blood
diggdt(1) = 0; %maternal blood
diggdt(2) = 0; %maternal blood
diggdt(3) = 0; %maternal blood
diggdt(4) = 0; %maternal blood

dy(5) = (params.k_up*y(1) - params.k_deg*y(5) + params.koff1*y(9) - params.kon1*y(5)*y(13))/params.v_stb; %STB endosomes unbound IgG1
dy(6) = (params.k_up*y(2) - params.k_deg*y(6) + params.koff2*y(10) - params.kon2*y(6)*y(13))/params.v_stb; %STB endosomes
dy(7) = (params.k_up*y(3) - params.k_deg*y(7) + params.koff3*y(11) - params.kon3*y(7)*y(13))/params.v_stb; %STB endosomes
dy(8) = (params.k_up*y(4) - params.k_deg*y(8) + params.koff4*y(12) - params.kon4*y(8)*y(13))/params.v_stb; %STB endosomes

dy(9) = (-params.koff1*y(9) + params.kon1*y(5)*y(13) - params.k_t*y(9))/params.v_stb; %concentration of params.params.fcRn-IgG1 compley
dy(10) = (-params.koff2*y(10) + params.kon2*y(6)*y(13)  - params.k_t*y(10))/params.v_stb; %concentration of params.params.fcRn-IgG2 compley
dy(11) = (-params.koff3*y(11) + params.kon3*y(7)*y(13)  - params.k_t*y(11))/params.v_stb; %concentration of params.params.fcRn-IgG3 compley
dy(12) = (-params.koff4*y(12) + params.kon4*y(8)*y(13)  - params.k_t*y(12))/params.v_stb; %concentration of params.params.fcRn-IgG4 compley

dy(13) =(-params.kon1*y(13)*y(5) + params.koff1*y(9) -params.kon2*y(13)*y(6) + params.koff2*y(10) ...
     - params.kon3*y(13)*y(7) + params.koff3*y(11) - params.kon4*y(13)*y(8) + params.koff4*y(12) ...
     + params.k_t*y(9) + params.k_t*y(10) + params.k_t*y(11) + params.k_t*y(12) + ...
     (2*params.fcrn_curve.p1*t + params.fcrn_curve.p2))/params.v_stb; %unbound params.params.fcRn in endosome

dy(14) = (params.k_t*y(9) + params.koff1b*y(18) - params.kon1b*y(14)*y(22))/params.v_str; %stroma unbound IgG1
dy(15) = (params.k_t*y(10) + params.koff2b*y(19) - params.kon2b*y(15)*y(22))/params.v_str; %stroma
dy(16) = (params.k_t*y(11) + params.koff3b*y(20) - params.kon3b*y(16)*y(22))/params.v_str; %stroma
dy(17) = (params.k_t*y(12) + params.koff4b*y(21) - params.kon4b*y(17)*y(22))/params.v_str; %stroma

dy(18) = (-params.koff1b*y(18) + params.kon1b*y(14)*y(22) - params.k_up*y(18))/params.v_str; %concentration of params.params.fcgRIIb-IgG1 compley
dy(19) = (-params.koff2b*y(19) + params.kon2b*y(15)*y(22)  - params.k_up*y(19))/params.v_str; %concentration of params.params.fcgRIIb-IgG2 compley
dy(20) = (-params.koff3b*y(20) + params.kon3b*y(16)*y(22)  - params.k_up*y(20))/params.v_str; %concentration of params.params.fcgRIIb-IgG3 compley
dy(21) = (-params.koff4b*y(21) + params.kon4b*y(17)*y(22)  - params.k_up*y(21))/params.v_str; %concentration of params.params.fcgRIIb-IgG4 compley

dy(22) =(-params.kon1b*y(22)*y(14) + params.koff1b*y(18) -params.kon2b*y(22)*y(15) + params.koff2b*y(19) ...
     - params.kon3b*y(22)*y(16) + params.koff3b*y(20) - params.kon4b*y(22)*y(17) + params.koff4b*y(21) ...
     + params.k_up*y(18) + params.k_up*y(19) + params.k_up*y(20) + params.k_up*y(21) + ...
     (2*params.fcgr2b_curve.p1*t + params.fcgr2b_curve.p2))/params.v_str; %unbound params.params.fcgRIIb on EC surface

dy(23) = (params.k_up*y(18) - params.k_t*y(23))/params.v_ec; %endothelial cells
dy(24) = (params.k_up*y(19) - params.k_t*y(24))/params.v_ec; %endothelial cells
dy(25) = (params.k_up*y(20) - params.k_t*y(25))/params.v_ec; %endothelial cells
dy(26) = (params.k_up*y(21) - params.k_t*y(26))/params.v_ec; %endothelial cells

dy(27) = (params.k_t*y(23) - params.d_ab*y(27))/params.v_f; %fetal blood
dy(28) = (params.k_t*y(24) - params.d_ab*y(28))/params.v_f; %fetal blood
dy(29) = (params.k_t*y(25) - params.d_ab*y(29))/params.v_f; %fetal blood
dy(30) = (params.k_t*y(26) - params.d_ab*y(30))/params.v_f; %fetal blood

end
