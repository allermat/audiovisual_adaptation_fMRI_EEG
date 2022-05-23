function [roi, AIC] = group_NF_FFX(dataPath)

% Psychometric function
PF = @PAL_CumulativeNormal;

% Stimulus locations
StimLevels = [-12 -5 -2 0 2 5 12];

% Load data
load(dataPath); 

% Number of regions
nROIs = numel(out.ROI);

% Set seed for reproducibility
rng('default');

% Format data
for i=1:nROIs
    for ss=1:numel(out.data)
        for j=1:3
            OutOfNum{i}(j,:,ss) = out.data(ss).OutOfNum{i,j};
            NumPos{i}(j,:,ss) = out.data(ss).NumPos{i,j};
        end
    end
    OutOfNum{i} = sum(OutOfNum{i}, 3);
    NumPos{i} = sum(NumPos{i}, 3);
end

for i=1:nROIs
    % Initial fit to get parameter estimates
    for j=1:3
        paramsInit(j,:) = fit_PF(PF, StimLevels, NumPos{i}(j,:), ...
            OutOfNum{i}(j,:), [], get_PAL_opt_precise);
    end
    
    % Eta parameter for betabinomial model
    eta = 0.15;
        
    %----- Model comparison between recalibration and static models
    
    % Fit lesser model (static model)
    thetasInit = [mean(paramsInit(:,1)) mean(paramsInit(:,2)) paramsInit(1,4)];
    paramsIDmatrix = [1 2 3 3; 1 2 3 3; 1 2 3 3];
    roi(i).mdl_lesser = multifit_PF(PF, 'lesser model', thetasInit, ...
        paramsIDmatrix, StimLevels, NumPos{i}, OutOfNum{i}, [0 0.45], eta);
    
    % Fit fuller model (recalibration model)
    thetasInit = [paramsInit(:,1)' mean(paramsInit(:,2)) paramsInit(1,4)];
    paramsIDmatrix = [1 4 5 5; 2 4 5 5; 3 4 5 5];
    roi(i).mdl_fuller = multifit_PF(PF, 'fuller model', thetasInit, ...
        paramsIDmatrix, StimLevels, NumPos{i}, OutOfNum{i}, [0 0.45], eta);
    
    % Calculate AIC (Akaike's Information Criterion) with adjustment 
    % suggested by Penny et al 2004
    AIC(i,1) = roi(i).mdl_fuller.LL - 6; % 3 PSE, 1 slope, 1 gamma/lapse and 1 eta parameter
    AIC(i,2) = roi(i).mdl_lesser.LL - 4; % 1 PSE, 1 slope, 1 gamma/lapse and 1 eta parameter
    
    %----- Calculate goodness of fit
        
    nBootstraps = 5000;
    for b=1:nBootstraps
        fprintf('Bootstrap %d...\n', b)
        
        % Generate simulated data
        for j=1:3
            NumPos_Sim(j,:) = DM_PF_SimulateObserverParametric_Beta(...
                roi(i).mdl_fuller.params(j,:), StimLevels, OutOfNum{i}(j,:), ...
                PF, roi(i).mdl_fuller.eta);
        end
        
        % Fit bootstrapped model
        roi(i).mdl_bootstrap(b) = multifit_PF(PF, sprintf('bootstrap %d', b), ...
            thetasInit, paramsIDmatrix, StimLevels, NumPos_Sim, OutOfNum{i}, ...
            [0 0.45], eta);
    end
    
    % Calculate p-value
    roi(i).mdl_fuller.pDev = (1+sum(cat(1, roi(i).mdl_bootstrap.Dev) > ...
        roi(i).mdl_fuller.Dev)) / (nBootstraps+1);
end