function [out, exp_r, xp, pxp, AIC] = group_PF_RFX(modality) 

% Psychometric function
PF = @PAL_CumulativeNormal;

% Stimulus locations
StimLevels = [-12 -5 -2 0 2 5 12];

% Load data
load([modality '_behav_all.mat']);

% Number of subjects
nsubjects = numel(out);

% Set seed for reproducibility
rng('default');

for ss=1:nsubjects
    fprintf('Subject %d...\n', ss)
    
    % Pretest
    NumPos{ss}(1,:) = sum(out(ss).pre.NumPos, 1);
    OutOfNum{ss}(1,:) = sum(out(ss).pre.OutOfNum, 1);
    
    % Posttest
    if ss == 1 && strcmp(modality, 'EEG')
        NumPos{ss}(2,:) = out(ss).post.NumPos(1,:); % radapt
        OutOfNum{ss}(2,:) = out(ss).post.OutOfNum(1,:);
        NumPos{ss}(3,:) = out(ss).post.NumPos(2,:); % ladapt
        OutOfNum{ss}(3,:) = out(ss).post.OutOfNum(2,:);
    else
        NumPos{ss}(2,:) = sum(out(ss).post.NumPos(1:2,:), 1); % radapt
        OutOfNum{ss}(2,:) = sum(out(ss).post.OutOfNum(1:2,:), 1);
        NumPos{ss}(3,:) = sum(out(ss).post.NumPos(3:4,:), 1); % ladapt
        OutOfNum{ss}(3,:) = sum(out(ss).post.OutOfNum(3:4,:), 1);
    end
    
    % Initial fit to get parameter estimates
    for i=1:3
        paramsInit(i,:) = fit_PF(PF, StimLevels, NumPos{ss}(i,:), ...
            OutOfNum{ss}(i,:), [], get_PAL_opt_precise);
    end
    
    % Eta parameter for betabinomial model
    eta = 0.15;
    
    %----- Model comparison between recalibration and static models
    
    % Fit lesser model (static model)
    thetasInit = [mean(paramsInit(:,1)) mean(paramsInit(:,2)) paramsInit(1,4)];
    paramsIDmatrix = [1 2 3 3; 1 2 3 3; 1 2 3 3];
    out(ss).mdl_lesser = multifit_PF(PF, 'lesser model', thetasInit, ...
        paramsIDmatrix, StimLevels, NumPos{ss}, OutOfNum{ss}, [0 0.1], eta);
    
    % Fit fuller model (recalibration model)
    thetasInit = [paramsInit(:,1)' mean(paramsInit(:,2)) paramsInit(1,4)];
    paramsIDmatrix = [1 4 5 5; 2 4 5 5; 3 4 5 5];
    out(ss).mdl_fuller = multifit_PF(PF, 'fuller model', thetasInit, ...
        paramsIDmatrix, StimLevels, NumPos{ss}, OutOfNum{ss}, [0 0.1], eta);

    % Calculate AIC (Akaike's Information Criterion) with adjustment 
    % suggested by Penny et al 2004
    AIC(ss,1) = out(ss).mdl_fuller.LL - 6; % 3 PSE, 1 slope, 1 gamma/lapse and 1 eta parameter
    AIC(ss,2) = out(ss).mdl_lesser.LL - 4; % 1 PSE, 1 slope, 1 gamma/lapse and 1 eta parameter
    
    %----- Calculate goodness of fit
        
    nBootstraps = 5000;
    for b=1:nBootstraps
        fprintf('Bootstrap %d...\n', b)
        
        % Generate simulated data
        for i=1:3
            NumPos_Sim(i,:) = DM_PF_SimulateObserverParametric_Beta(...
                out(ss).mdl_fuller.params(i,:), StimLevels, OutOfNum{ss}(i,:), ...
                PF, out(ss).mdl_fuller.eta);
        end
        
        % Fit bootstrapped model
        out(ss).mdl_bootstrap(b) = multifit_PF(PF, sprintf('bootstrap %d', b), ...
            thetasInit, paramsIDmatrix, StimLevels, NumPos_Sim, OutOfNum{ss}, ...
            [0 0.1], eta);
    end
    
    % Calculate p-value
    out(ss).mdl_fuller.pDev = (1+sum(cat(1, out(ss).mdl_bootstrap.Dev) > ...
        out(ss).mdl_fuller.Dev)) / (nBootstraps+1);
end

% Calculate protected exceedance probability
[alpha, exp_r, xp, pxp] = spm_BMS(AIC);