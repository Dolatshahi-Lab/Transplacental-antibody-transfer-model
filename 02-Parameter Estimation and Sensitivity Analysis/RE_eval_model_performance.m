function [passedtotalModelRuns, passedtotalParams, failedtotalModelRuns, failedtotalParams] = RE_eval_model_performance(calibrationInput, modelRuns, passCriteria)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
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
% Evaluates model performance by comparing model run outcomes against
% experimental data provided by the user.
%
% INPUTS
%   - calibrationInput: structure that holds user input to calibration
%   program.
%   - modelRuns: structure that holds modelRuns output from the lhs_ode_run
%   function as well as other information such as paramMatrix.
%
% OUTPUTS
%   - passedtotalModelRuns -- a matrix that holds the runs that passed
%   following comparison to user-provided data.
%   - passedtotalParams:  a matrix holding the parameter values for the
%   runs that passed following comparison to user-provided data.
%   - failedtotalModelRuns:  a matrix that holds the runs that failed
%   following comparison to user-provided data. 
%   - failedtotalParams: a matrix that holds the parameter values for the
%   runs that failed following comparison to user-provided data. 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%name some tempvars for model Output:
% preyOutput = modelRuns.odeOut{1,1}*1000; % timesteps x number of runs for output 1 -- this is multipled by 1000 because the experimental data is in the thousands
% predOutput = modelRuns.odeOut{2,1}*1000; % timesteps x number of runs for output 2 -- this is multipled by 1000 because the experimental data is in the thousands
igg1Output = modelRuns.odeOut{27,1}; % timesteps x number of runs for output 1
igg2Output = modelRuns.odeOut{28,1}; % timesteps x number of runs for output 2
igg3Output = modelRuns.odeOut{29,1}; % timesteps x number of runs for output 1
igg4Output = modelRuns.odeOut{30,1}; % timesteps x number of runs for output 2

%name some tempvars for experimental data:
% hareDS = calibrationInput.HareDS;
% lynxDS = calibrationInput.LynxDS;
igg1DS = calibrationInput.igg1DS;
igg2DS = calibrationInput.igg2DS;
igg3DS = calibrationInput.igg3DS;
igg4DS = calibrationInput.igg4DS;

% how does the output perform against dataset?
% [passedPreyRuns, runIndxPreyPassed, runIndxPreyFailed] = find_passed(preyOutput, hareDS,passCriteria);
% [passedPredRuns, runIndxPredPassed, runIndxPredFailed] = find_passed(predOutput, lynxDS,passCriteria);
[~, runIndxIgG1Passed, runIndxIgG1Failed] = find_passed(igg1Output(modelRuns.odeSettings.time_index,:), igg1DS,passCriteria);
[~, runIndxIgG2Passed, runIndxIgG2Failed] = find_passed(igg2Output(modelRuns.odeSettings.time_index,:), igg2DS,passCriteria);
[~, runIndxIgG3Passed, runIndxIgG3Failed] = find_passed(igg3Output(modelRuns.odeSettings.time_index,:), igg3DS,passCriteria);
[~, runIndxIgG4Passed, runIndxIgG4Failed] = find_passed(igg4Output(modelRuns.odeSettings.time_index,:), igg4DS,passCriteria);

%find the columns that are in common between the two conditions (prey data
%and predator data)
CommonIndxPassed = intersect(runIndxIgG1Passed,runIndxIgG2Passed);
CommonIndxPassed = intersect(CommonIndxPassed,runIndxIgG3Passed);
CommonIndxPassed = intersect(CommonIndxPassed,runIndxIgG4Passed);

CommonIndxFailed = union(runIndxIgG1Failed, runIndxIgG2Failed);
CommonIndxFailed = union(CommonIndxFailed, runIndxIgG3Failed);
CommonIndxFailed = union(CommonIndxFailed, runIndxIgG4Failed);

%subset runs, LHS matrix, and param matrix to include only the ones that
%passed:
passedtotalParams = modelRuns.paramMatrix(CommonIndxPassed,:);
passedtotalModelRuns(:,:,1) = igg1Output(:,CommonIndxPassed);
passedtotalModelRuns(:,:,2) = igg2Output(:,CommonIndxPassed);
passedtotalModelRuns(:,:,3) = igg3Output(:,CommonIndxPassed);
passedtotalModelRuns(:,:,4) = igg4Output(:,CommonIndxPassed);

%subset runs and param matrix to also include those that did not pass:
failedtotalParams = modelRuns.paramMatrix(CommonIndxFailed,:);
failedtotalModelRuns(:,:,1) = igg1Output(:,CommonIndxFailed);
failedtotalModelRuns(:,:,2) = igg2Output(:,CommonIndxFailed);
failedtotalModelRuns(:,:,3) = igg3Output(:,CommonIndxFailed);
failedtotalModelRuns(:,:,4) = igg4Output(:,CommonIndxFailed);

end


function [passedModelRuns, RunIndxThatPassed, RunIndxThatFailed] = find_passed(modelOutput, expData, passCriteria)

% Customized pass criteria - simulation can't be higher than the upper
% limit at time point 1 or final, and can't be lower than the lower limit at final
% time point
LogicalArrayCondition1 = modelOutput(end,:) >= expData(end,1) / passCriteria;
LogicalArrayCondition2 = modelOutput([1,end],:) <= expData([1,end],2) * passCriteria;
LogicalArraytot = LogicalArrayCondition1 & LogicalArrayCondition2;

% cols with all ones --  meaning they passed the criteria
RunIndxThatPassed = find(all(LogicalArraytot==1)); 

%cols with all zeroes -- meaning they failed the criteria:
RunIndxThatFailed = find(~all(LogicalArraytot));

% %logical indexing to quickly check if the runs pass through the window of
% %our data (time or divided by our passing criteria). 
% LogicalArrayCondition1 = modelOutput(:,:) >= expData(:,1) / passCriteria;
% LogicalArrayCondition2 = modelOutput(:,:) <= expData(:,2) * passCriteria;
% LogicalArraytot = LogicalArrayCondition1 & LogicalArrayCondition2;

% % cols with all ones --  meaning they passed the criteria
% RunIndxThatPassed = find(all(LogicalArraytot==1)); 
% 
% %cols with all zeroes -- meaning they failed the criteria:
% RunIndxThatFailed = find(~all(LogicalArraytot));

%use RunIndxThatPassed to pull those columns from the total runs matrix:
passedModelRuns = modelOutput(:,RunIndxThatPassed);

end
