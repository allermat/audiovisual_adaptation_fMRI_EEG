% Figure 5
clearvars;
close all;
% Panel A: BOLD responses
dataPath = fullfile(get_path('project'),'results','data',...
                    'bold_fMRI.mat');
if ~exist(dataPath,'file')
    out = run_bold_fMRI('hemisphere','mean_lhrh');
    save(dataPath,'out');
end
bold_table = figure_BOLD(dataPath);

% Panel B: linear mixed effects (LME) models
dataPath = fullfile(get_path('project'),'results','data',...
                    'lme_fMRI_pre_lh_and_rh.mat');
if ~exist(dataPath,'file')
    out = run_lme_fMRI_pre();
    save(dataPath,'out');
end
logBF_table_lme_pre = figure_lme_fMRI_pre(dataPath);

dataPath = fullfile(get_path('project'),'results','data',...
                    'lme_fMRI_prepost_lh_and_rh.mat');
if ~exist(dataPath,'file')
    out = run_lme_fMRI_prepost();
    save(dataPath,'out');
end
logBF_table_lme_prepost = figure_lme_fMRI_prepost(dataPath);

% Panel C: pattern component modelling (PCM)
dataPath = fullfile(get_path('project'),'results','data',...
                    'pcm_fMRI_pre_hemi_dec_fixedEff_noEucnorm.mat');
if ~exist(dataPath,'file')
    [out,M] = run_pcm_fMRI_pre('runEffect','fixed');
    save(dataPath,'out','M');
end
logBF_table_pcm_pre = figure_PCM_fMRI_pre(dataPath);

dataPath = fullfile(get_path('project'),'results','data',...
                    'pcm_fMRI_prepost_hemi_dec_randomEff_noEucnorm.mat');
if ~exist(dataPath,'file')
    [out,M] = run_pcm_fMRI_prepost('runEffect','random','iseucnorm',false);
    save(dataPath,'out','M');
end
logBF_table_pcm_prepost = figure_PCM_fMRI_prepost(dataPath);