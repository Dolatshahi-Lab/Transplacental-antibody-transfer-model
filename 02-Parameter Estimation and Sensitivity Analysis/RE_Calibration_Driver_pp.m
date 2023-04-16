function RE_Calibration_Driver_pp(lhsODESettingsFileName, seedGenAndNum)
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% This script sourced from:  http://malthus.micro.med.umich.edu/CaliPro/
% Joslyn, L. et al. CaliPro:  A Calibration Protocol That Utilizes Parameter
% Density Estimation to Explore Parameter Space and Calibrate Complex 
% Biological Models.  Cellular and Molecular Bioengineering (2021).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% And adapted by:  Erdogan, R. A quantitative mechanistic model reveals
% key determinants of maternal-fetal IgG transfer with implications for 
% prenatal immunization (2023).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is a Driver file for an automated, model-agnostic calibration
% strategy of complex systems that require unconventional objective
% function evaluations (i.e., non-continuous objective functions) and
% utilitzes an LHS sampling scheme for sweeping global parameter space.
% However, this process could be modified for any sampling scheme, and is
% not dependent on LHS sampling.  For example, Sobol sequencing could be
% used.
%
% CaliPro (as implemented below) is currently used with an ODE model that
% utilizes LHS sampling.  The sampling and model specification setup
% follows the process our lab has used for several years now, but the
% general approach of CaliPro should work for any model, including
% stochastic models like ABMs, and we hope to generalize this framework for
% future use.
%
% The steps of the calibration method are the following:
%       1. Provide the Program with the Necessary Information to Run.
%            - initially, this includes the model itself, the starting
%              parameter ranges, any data to compare against (experimental
%              data, for example) and the termination condition.
%       2. Perform Sampling According to the Sampling Scheme (LHS). 
%       3. Run the Model 
%              - steps 2 and 3 might be combined in using the
%              Predator-Prey model system evaluated using our lhs_ode lab
%              setup
%        4. Evaluate the Performance of the Model Against Experimental Data.
%            - Create Passed Model Runs Data structure
%            - Create Passed Model Runs Params Data structure
%            - if performace satisfies termination condition, then leave 
%              the calibration program.
%       5. Perform Density (Distribution) Comparisons between Passed Model
%          Parameter and Total Model Parameter Ranges.
%            - Evaluate the densest portions
%              of the passed parameter set ranges. 
%            - The range of these densest regions are the new parameter
%              ranges for the next iteration of the calibration program. 
%            - AlternateRun designates if the user wishes to use the
%            alternate method of comparing passed runs vs failed runs via
%            substitution of the latter from the former.  This substitution
%            gives a range of parameter values wherein the runs mostly
%            passed. 
% 
% Repeat this process (particularly steps 2-5) until condition is
% satisfied. As an arbitrary example, the condition might be satisfied when
% 90% of individiuals have an infection where all granulomas can control
% the bacterial growth. 
% 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
%     1. Provide the Program with the Necessary Information to Run.
% -------------------------------------------------------------------------

% model information is stored in the model settings file 
% Load the settings file - everything returned in the odeSettings structure.
odeSettings = lhs_ode_get_run_settings_new(lhsODESettingsFileName);

% if seed information was input, then assign that information to the seed
% generator and seed number.
if exist('seedGenAndNum', 'var')
    rng(seedGenAndNum);
end
% save the current generator settings for reproduciblity.
% this approach allows for long-term reproducibility by saving not only the
% seed, but also the seed generator.
odeSettings.seedGenSettings = rng;

% eventually, load the calibration settings file - everything returned in
% the calibration structure. For now, just put the information inside the
% driver file

calibrationInput.HDR = false;

%How long should this run for unsuccesfully?
calibrationInput.maxIter = 20;

%What is the percentage of runs that need to satisfy our criteria for the
%termination condition to be executed?
calibrationInput.terminate = 0.95;

%What is the portion of parameters that we should be selecting
%within the parameter ranges at the end of each iteration? 
calibrationInput.portionParams = 0.85; 

%how flexible is the passing criteria? 
calibrationInput.passCriteria = 5;

%what are the names of the experimental datasets?  Datasets to compare
%against model outcomes  
calibrationInput.igg1DS = csvread('igg1data.csv');
calibrationInput.igg2DS = csvread('igg2data.csv');
calibrationInput.igg3DS = csvread('igg3data.csv');
calibrationInput.igg4DS = csvread('igg4data.csv');

calibrationInput.passCriteriaInitial = 5;
calibrationInput.passCriteriaNonInitial = 1.25;

%set up the array that holds all of the calibration results for each run:
CalibrationOutputArray = {};
%loop through, until condition is satisfied:
runIteration = 1;
terminateCalib = false; 
while ~terminateCalib
    %for runIteration = 1:1
    % -------------------------------------------------------------------------
    %     2. Perform Sampling According to the Sampling Scheme (LHS). 
    % -------------------------------------------------------------------------
    % -------------------------------------------------------------------------
    %     3. Run the Model.
    % -------------------------------------------------------------------------
    %for this model framework, LHS sampling and model execution happen through
    %calling the lhs_ode_run_new_c function.  the lhs_ode_run_new function has
    %been altered to work better within the calibration framework.
    modelRuns = RE_lhs_ode_run_new_c(odeSettings,runIteration);

    % -------------------------------------------------------------------------
    %     4. Evaluate the Performance of the Model Against Experimental Data.
    % -------------------------------------------------------------------------

    %make fit criteria stricter after initial run:
    passCriteria = calibrationInput.passCriteriaInitial;
    if runIteration > 2
        passCriteria = calibrationInput.passCriteriaNonInitial;
    end

    % Compare model output against the experimental data
    [passedtotalModelRuns, passedtotalParams, failedtotalModelRuns, failedtotalParams] = RE_eval_model_performance(calibrationInput, modelRuns, passCriteria);
    % save output/plot results
    modelRuns.passedtotalModelRuns = passedtotalModelRuns;
    modelRuns.passedtotalParams = passedtotalParams;
    modelRuns.failedtotalModelRuns = failedtotalModelRuns;
    modelRuns.failedtotalParams = failedtotalParams;

    % evaluate the performance according to termination conditions
    %leave loop once condition is satisfied:
    if (size(passedtotalModelRuns,2) / odeSettings.NR) > calibrationInput.terminate
        disp('Calibration Complete According to User Specifications')
        CalibrationOutputArray{runIteration} = modelRuns;
        terminateCalib = true;
        break
    end 

    % -------------------------------------------------------------------------
    %     5. Perform Density (Distribution) Comparisons between Passed Model
    %        Parameter and Total Model Parameter Ranges.
    % -------------------------------------------------------------------------   
    %   if user specfied,             
    if ~(calibrationInput.HDR) %find the param ranges using alternate method via subtraction method
        [hdr_param1, density1, xi1, xdensityFail1, xFail1, xdensitySub1] = find_ads(passedtotalParams(:,1), failedtotalParams(:,1));
        [hdr_param2, density2, xi2, xdensityFail2, xFail2, xdensitySub2] = find_ads(passedtotalParams(:,2), failedtotalParams(:,2));
        [hdr_param3, density3, xi3, xdensityFail3, xFail3, xdensitySub3] = find_ads(passedtotalParams(:,3), failedtotalParams(:,3));
        [hdr_param4, density4, xi4, xdensityFail4, xFail4, xdensitySub4] = find_ads(passedtotalParams(:,4), failedtotalParams(:,4));
        [hdr_param5, density5, xi5, xdensityFail5, xFail5, xdensitySub5] = find_ads(passedtotalParams(:,5), failedtotalParams(:,5));
        [hdr_param6, density6, xi6, xdensityFail6, xFail6, xdensitySub6] = find_ads(passedtotalParams(:,6), failedtotalParams(:,6));
        [hdr_param7, density7, xi7, xdensityFail7, xFail7, xdensitySub7] = find_ads(passedtotalParams(:,7), failedtotalParams(:,7));
        [hdr_param8, density8, xi8, xdensityFail8, xFail8, xdensitySub8] = find_ads(passedtotalParams(:,8), failedtotalParams(:,8));
        [hdr_param9, density9, xi9, xdensityFail9, xFail9, xdensitySub9] = find_ads(passedtotalParams(:,9), failedtotalParams(:,9));
        [hdr_param10, density10, xi10, xdensityFail10, xFail10, xdensitySub10] = find_ads(passedtotalParams(:,10), failedtotalParams(:,10));

    else %find hdr for each param in passedtotal params. Return failed info too   
        [hdr_param1, density1, xi1, xdensityFail1, xFail1] = find_hdr(passedtotalParams(:,1), failedtotalParams(:,1), calibrationInput.portionParams);
        [hdr_param2, density2, xi2, xdensityFail2, xFail2] = find_hdr(passedtotalParams(:,2), failedtotalParams(:,2), calibrationInput.portionParams);
        [hdr_param3, density3, xi3, xdensityFail3, xFail3] = find_hdr(passedtotalParams(:,3), failedtotalParams(:,3), calibrationInput.portionParams);
        [hdr_param4, density4, xi4, xdensityFail4, xFail4] = find_hdr(passedtotalParams(:,4), failedtotalParams(:,4), calibrationInput.portionParams);
        [hdr_param5, density5, xi5, xdensityFail5, xFail5] = find_hdr(passedtotalParams(:,5), failedtotalParams(:,5), calibrationInput.portionParams);
        [hdr_param6, density6, xi6, xdensityFail6, xFail6] = find_hdr(passedtotalParams(:,6), failedtotalParams(:,6), calibrationInput.portionParams);
        [hdr_param7, density7, xi7, xdensityFail7, xFail7] = find_hdr(passedtotalParams(:,7), failedtotalParams(:,7), calibrationInput.portionParams);
        [hdr_param8, density8, xi8, xdensityFail8, xFail8] = find_hdr(passedtotalParams(:,8), failedtotalParams(:,8), calibrationInput.portionParams);
        [hdr_param9, density9, xi9, xdensityFail9, xFail9] = find_hdr(passedtotalParams(:,9), failedtotalParams(:,9), calibrationInput.portionParams);
        [hdr_param10, density10, xi10, xdensityFail10, xFail10] = find_hdr(passedtotalParams(:,10), failedtotalParams(:,10), calibrationInput.portionParams);
    end

    %load the new parameter ranges, based on the output of the find_hdr or
    %find_hdr_alt.  If the number of parameters is larger, there could be some
    %sort of function or (at the very least) a for loop that does this for all
    %the parameters. Could also include a for loop for the if else statement
    %above.
    odeSettings.parameters = ...
    {
    { 'k_{up}', 'u', hdr_param1(1), hdr_param1(2)} ... 0.5, 0.7 } ...%
    { 'k_{t}',  'u', hdr_param2(1), hdr_param2(2)} ...'u', 0.02, 0.035  } ...%
    { 'FcRn', 'u', hdr_param3(1), hdr_param3(2)} ...0.6, 0.9  } ...%
    { 'FcgRIIb', 'u', hdr_param4(1), hdr_param4(2)} ...'u', 0.02, 0.03 } ...%
    { 'v_M', 'u', hdr_param5(1), hdr_param5(2)} ...'u', 0.02, 0.03 } ...%
    { 'v_F', 'u', hdr_param6(1), hdr_param6(2)} ...'u', 0.02, 0.03 } ...%
    { 'v_{STR}', 'u', hdr_param7(1), hdr_param7(2)} ...'u', 0.02, 0.03 } ...%
    { 'v_{STB}', 'u', hdr_param8(1), hdr_param8(2)} ...'u', 0.02, 0.03 } ...%
    { 'v_{EC}', 'u', hdr_param9(1), hdr_param9(2)} ...'u', 0.02, 0.03 } ...%
    { 'k_{deg}', 'u', hdr_param10(1), hdr_param10(2)} ...'u', 0.02, 0.03 } ...%
    };

    %update modelRuns structure with the Density values for those parameters
    %whose simulation resulted in a pass/fail, the parameter range for those
    %runs and, if specified, the density values that represent the subtraction
    %of passed parameters from failed parameters
    modelRuns.ParamsDensityPass = [density1; density2; density3; density4; density5; density6; density7; density8; density9; density10];
    modelRuns.ParamsDensityRangePass = [xi1; xi2; xi3; xi4; xi5; xi6; xi7; xi8; xi9; xi10];
    modelRuns.ParamsDensityFail = [xdensityFail1; xdensityFail2; xdensityFail3; xdensityFail4;  xdensityFail5;  xdensityFail6;  xdensityFail7;  xdensityFail8;  xdensityFail9;  xdensityFail10];
    modelRuns.ParamsDensityRangeFail = [xFail1; xFail2; xFail3; xFail4;  xFail5;  xFail6;  xFail7;  xFail8;  xFail9;  xFail10];
    if ~(calibrationInput.HDR) 
        modelRuns.DensitySub = [xdensitySub1; xdensitySub2; xdensitySub3; xdensitySub4;  xdensitySub5;  xdensitySub6;  xdensitySub7;  xdensitySub8;  xdensitySub9;  xdensitySub10];
    end
    modelRuns.hdrRange = [hdr_param1; hdr_param2; hdr_param3; hdr_param4; hdr_param5; hdr_param6; hdr_param7; hdr_param8; hdr_param9; hdr_param10];

    %update CalibrationOutputArray with all the information from this run
    %iteration in the calibration program. 
    CalibrationOutputArray{runIteration} = modelRuns;

    %update the runIteration counter
    runIteration = runIteration + 1;

    % evaluate the performance according to termination conditions
    %leave loop once condition is satisfied:
    if (runIteration > calibrationInput.maxIter)
        disp('Calibration Met Max Iteration Count. Evaluate initial parameter ranges or model construction')
        terminateCalib = 1;
    end 

end  %end while loop that performs the calibration

% Assign variables to the base workspace so they will be saved after
% leaving the calibration.
assignin('base', 'CalibrationOutputArray', CalibrationOutputArray);
assignin('base', 'calibrationInput', calibrationInput);
assignin('base', 'passedtotalModelRuns', passedtotalModelRuns); %Not necessary long-term, helpful for debug right now tho
assignin('base', 'passedtotalParams', passedtotalParams); %Not necessary long-term, helpful for debug right now tho
assignin('base', 'odeSettings', odeSettings); %Not necessary long-term, helpful for debug right now tho
end %end calibration function
