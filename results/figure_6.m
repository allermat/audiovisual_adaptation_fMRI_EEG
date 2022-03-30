% Figure 6
clearvars;
close all;
% Panel A, B: ERPs, decoding accuracy and recalibration index
figure_EEG_results;

% Panel C: EEG PCM
dataPath = fullfile(get_path('project'),'results','data',...
                    'pcm_EEG_pre_hemi_dec_randomEff_noEucnorm.mat');
if ~exist(dataPath,'file')
    [out,M] = run_pcm_EEG_pre('runEffect','random');
    save(dataPath,'out','M');
end
logBF_eeg_pcm_pre = figure_PCM_EEG_pre(dataPath);

dataPath = fullfile(get_path('project'),'results','data',...
                    'pcm_EEG_prepost_hemi_dec_randomEff_noEucnorm.mat');
if ~exist(dataPath,'file')
    [out,M] = run_pcm_EEG_prepost('runEffect','random');
    save(dataPath,'out','M');
end
logBF_eeg_pcm_prepost = figure_PCM_EEG_prepost(dataPath);

% Panel D: EEG PCM with fMRI ROIs as predictor models
dataPath = fullfile(get_path('project'),'results','data',...
                    'pcm_fMRI_to_EEG_prepost_randomEff_noEucnorm.mat');
if ~exist(dataPath,'file')
    [out,M] = run_pcm_fMRI_to_EEG_prepost('runEffect','random');
    save(dataPath,'out','M');
end
logBF_eeg_fmri_pcm = figure_PCM_fMRI_to_EEG_prepost(dataPath);