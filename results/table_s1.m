% Table S1
clearvars;
close all;
% Psychophysics
dataPath = fullfile(get_path('project'),'results','data',...
    'PP_behav_all.mat');
if ~exist(dataPath,'file')
    out = run_behav_all('psychophysics');
    save(dataPath,'out');
end
behav_table_pp = table_behav_all(dataPath);

% fMRI
dataPath = fullfile(get_path('project'),'results','data',...
    'fMRI_behav_all.mat');
if ~exist(dataPath,'file')
    out = run_behav_all('fMRI');
    save(dataPath,'out');
end
behav_table_fmri = table_behav_all(dataPath);

% EEG
dataPath = fullfile(get_path('project'),'results','data',...
    'EEG_behav_all.mat');
if ~exist(dataPath,'file')
    out = run_behav_all('EEG');
    save(dataPath,'out');
end
behav_table_eeg = table_behav_all(dataPath);