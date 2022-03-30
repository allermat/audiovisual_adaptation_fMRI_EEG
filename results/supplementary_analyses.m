%% Supplementary Analyses for Nature Communications
clearvars;
close all;
%% LME - Using only the pretest data
% With decisional model 
dataPath = fullfile(get_path('project'),'results','data',...
                    'lme_fMRI_pre_lh_and_rh.mat');
if ~exist(dataPath,'file')
    out = run_lme_fMRI_pre();
    save(dataPath,'out');
end
logBF_table_lme_pre = figure_lme_fMRI_pre(dataPath);

% Adding decisional choice to the winning model in each ROI to see if it
% increases evidence (The first model is always constant, as that is the
% null model as well). 
dataPath = fullfile(get_path('project'),'results','data',...
                    'lme_fMRI_pre_lh_and_rh_winning_pmf.mat');
if ~exist(dataPath,'file')
    modelNames = {{'Const','Pmf+Const'},...
              {'Const','Sp+Const','Pmf+Sp+Const'},...
              {'Const','Pmf+Const'},...
              {'Const','Dec+Const','Pmf+Dec+Const'},...
              {'Const','Dec+Const','Pmf+Dec+Const'}};
    out = run_lme_fMRI_pre_custom(modelNames,'modelColors',modelColors,...
                                  'hemisphere','lh_and_rh');
    save(dataPath,'out');
end

colors    = {[.7 0 0],...      % red
             [0 0 .7],...      % blue
             [.9 .6 0],...     % orange
             [0 0.6 0.6],...   % cyan
             [0.5 0 0.5],...   % purple
             [0.2 0.6 0.2]};   % green
modelColors = {{'k',colors{6}},...
               {colors{1},colors{6}},...
               {'k',colors{6}},...
               {colors{2},colors{6}},...
               {colors{2},colors{6}}};
logBF_table_lme_pre_wining_pmf = figure_lme_fMRI_custom(dataPath,...
                                    'modelColors',modelColors);

%% LME - on pre and posttest data
% With decisional model
dataPath = fullfile(get_path('project'),'results','data',...
                    'lme_fMRI_prepost_lh_and_rh.mat');
if ~exist(dataPath,'file')
    out = run_lme_fMRI_prepost();
    save(dataPath,'out');
end
logBF_table_lme_prepost = figure_lme_fMRI_prepost(dataPath);

% Adding decisional choice to the winning model in each ROI to see if it
% increases evidence (The first model is always constant, as that is the
% null model as well).
dataPath = fullfile(get_path('project'),'results','data',...
                    'lme_fMRI_prepost_lh_and_rh_winning_pmf.mat');
if ~exist(dataPath,'file')
    modelNames = {{'Const','PmfRecal+Const'},...
              {'Const','SpRecal+Dec+Const','PmfRecal+SpRecal+Dec+Const'},...
              {'Const','Sp+DecRecal+Const','PmfRecal+Sp+DecRecal+Const'},...
              {'Const','Sp+DecRecal+Const','PmfRecal+Sp+DecRecal+Const'},...
              {'Const','Sp+DecRecal+Const','PmfRecal+Sp+DecRecal+Const'}};
    out = run_lme_fMRI_prepost_custom(modelNames,'hemisphere','lh_and_rh');
    save(dataPath,'out');
end
colors    = {[.7 0 0],...      % red
             [0 0 .7],...      % blue
             [.9 .6 0],...     % orange
             [0 0.6 0.6],...   % cyan
             [0.5 0 0.5],...   % purple
             [0.2 0.6 0.2]};   % green
modelColors = {{'k',colors{6}},...
               {colors{2},colors{6}},...
               {colors{3},colors{6}},...
               {colors{3},colors{6}},...
               {colors{3},colors{6}}};
logBF_table_lme_prepost_wining_pmf = figure_lme_fMRI_custom(dataPath,...
                                         'modelColors',modelColors);

%% PCM - Using only the pretest data
dataPath = fullfile(get_path('project'),'results','data',...
                    'pcm_fMRI_pre_hemi_dec_fixedEff_noEucnorm.mat');
if ~exist(dataPath,'file')
    [out,M] = run_pcm_fMRI_pre('runEffect','fixed');
    save(dataPath,'out','M');
end
logBF_fMRI_pre = figure_PCM_fMRI_pre(dataPath);

% Taking the winning model in each ROI and adding the PMF model to them
idx_fixed = repmat(logical([1,1,1,0]'),5,1);

logBF_fMRI_pre.idx_fixed = idx_fixed;

temp = logBF_fMRI_pre(logBF_fMRI_pre.idx_fixed,:);
temp.input = cat(2,temp.model,num2cell(temp.mean));
winMdl_pre_dec_fixed = varfun(@(x) x(argmax(cell2mat(x(:,2))),:),temp,...
                          'InputVariables','input',...
                          'GroupingVariables',{'roi'});

% Taking the winning model in each ROI and adding the PMF model
dataPath = fullfile(get_path('project'),'results','data',...
                    'pcm_fMRI_pre_hemi_winning+pmf_fixedEff_noEucnorm.mat');
if ~exist(dataPath,'file')
    [out,M] = run_pcm_fMRI_pre_winning_pmf('runEffect','fixed');
    save(dataPath,'out','M');
end
logBF_fMRI_pre_winning_pmf = figure_PCM_fMRI_pre_winning_pmf(dataPath);

%% PCM - on pre- and posttest data
dataPath = fullfile(get_path('project'),'results','data',...
                    'pcm_fMRI_prepost_hemi_dec_randomEff_noEucnorm.mat');
if ~exist(dataPath,'file')
    [out,M] = run_pcm_fMRI_prepost('runEffect','random');
    save(dataPath,'out','M');
end
logBF_fMRI_prepost = figure_PCM_fMRI_prepost(dataPath);

% Getting highest model evidence separately for fixed and random runeffects 
idx_random = repmat(logical([1,1,1,1,0]'),5,1);
logBF_fMRI_prepost.idx_random = idx_random;
temp = logBF_fMRI_prepost(logBF_fMRI_prepost.idx_random,:);
temp.input = cat(2,temp.model,num2cell(temp.mean));
winMdl_prepost_dec_random = varfun(@(x) x(argmax(cell2mat(x(:,2))),:),temp,...
                           'InputVariables','input',...
                           'GroupingVariables',{'roi'});

% Taking the winning model in each ROI and adding the PMF model 
dataPath = fullfile(get_path('project'),'results','data',...
                    'pcm_fMRI_prepost_hemi_winning+pmf_randomEff_noEucnorm.mat');
if ~exist(dataPath,'file')
    [out,M] = run_pcm_fMRI_prepost_winning_pmf('runEffect','random');
    save(dataPath,'out','M');
end
logBF_fMRI_prepost_winning_pmf = figure_PCM_fMRI_prepost_winning_pmf(dataPath);

%% PCM EEG - 0n pretest data
dataPath = fullfile(get_path('project'),'results','data',...
                    'pcm_EEG_pre_hemi_dec_randomEff_noEucnorm.mat');
if ~exist(dataPath,'file')
    [out,M] = run_pcm_EEG_pre('runEffect','random');
    save(dataPath,'out','M');
end
logBF_EEG_pre = figure_PCM_EEG_pre(dataPath);

% Taking the winning model in each time window and adding the PMF model
temp = logBF_EEG_pre(~ismember(logBF_EEG_pre.model,'freechol'),:);
temp.input = cat(2,temp.model,num2cell(temp.mean));
winMdl_EEG_pre_dec = varfun(@(x) x(argmax(cell2mat(x(:,2))),:),temp,...
                          'InputVariables','input',...
                          'GroupingVariables',{'timewin'});

dataPath = fullfile(get_path('project'),'results','data',...
                    'pcm_EEG_pre_winning+pmf_randomEff_noEucnorm.mat');
if ~exist(dataPath,'file')
    [out,M] = run_pcm_EEG_pre_winning_pmf('runEffe','random');
    save(dataPath,'out','M');
end
logBF_EEG_pre_winning_pmf = figure_PCM_EEG_pre_winning_pmf(dataPath);

%% PCM EEG - 0n pre- and posttest data, revision improvements
dataPath = fullfile(get_path('project'),'results','data',...
                    'pcm_EEG_prepost_hemi_dec_randomEff_noEucnorm.mat');
if ~exist(dataPath,'file')
    [out,M] = run_pcm_EEG_prepost('runEffect','random');
    save(dataPath,'out','M');
end
logBF_EEG_prepost = figure_PCM_EEG_prepost(dataPath);

% Taking the winning model in each ROI and adding the PMF model to them
temp = logBF_EEG_prepost(~ismember(logBF_EEG_prepost.model,'freechol'),:);
temp.input = cat(2,temp.model,num2cell(temp.mean));
winMdl_EEG_prepost_dec = varfun(@(x) x(argmax(cell2mat(x(:,2))),:),temp,...
                          'InputVariables','input',...
                          'GroupingVariables',{'timewin'});

dataPath = fullfile(get_path('project'),'results','data',...
                    'pcm_EEG_prepost_winning+pmf_randomEff_noEucnorm.mat');
if ~exist(dataPath,'file')
    [out,M] = run_pcm_EEG_prepost_winning_pmf('runEffect','random',...
                                              'iseucnorm',false);
    save(dataPath,'out','M');
end
logBF_EEG_prepost_winning_pmf = figure_PCM_EEG_prepost_winning_pmf(dataPath);
