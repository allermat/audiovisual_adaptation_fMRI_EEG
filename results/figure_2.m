% Figure 2
clearvars;
close all;
% Panel A: Psychometric functions
dataPath = fullfile(get_path('project'),'results','data',...
    'PP_behav_all.mat');
if ~exist(dataPath,'file')
    out = run_behav_all('psychophysics');
    save(dataPath,'out');
end
pf_table_pp = figure_PF_prepost(dataPath);

dataPath = fullfile(get_path('project'),'results','data',...
    'fMRI_behav_all.mat');
if ~exist(dataPath,'file')
    out = run_behav_all('fMRI');
    save(dataPath,'out');
end
pf_table_fmri = figure_PF_prepost(dataPath);

dataPath = fullfile(get_path('project'),'results','data',...
    'EEG_behav_all.mat');
if ~exist(dataPath,'file')
    out = run_behav_all('EEG');
    save(dataPath,'out');
end
pf_table_eeg = figure_PF_prepost(dataPath);

% Panel B: PSE values
dataPath = fullfile(get_path('project'),'results','data',...
    'PP_behav_all.mat');
pse_table_prepost_pp = figure_PSE_prepost(dataPath);

dataPath = fullfile(get_path('project'),'results','data',...
    'fMRI_behav_all.mat');
pse_table_prepost_fmri = figure_PSE_prepost(dataPath);

dataPath = fullfile(get_path('project'),'results','data',...
    'EEG_behav_all.mat');
pse_table_prepost_eeg = figure_PSE_prepost(dataPath);

% Calculate bootstrap-based paired ttest on PSE values
data = pse_table_prepost_pp{2:3,4:end};
[p_uncorr,h_uncorr,~,~,obsStat] = mvpa.bootstrpOneSampleTtest(diff(data), ...
    10000, 'MCPsol', 'none')

data = pse_table_prepost_fmri{2:3,4:end};
[p_uncorr,h_uncorr,~,~,obsStat] = mvpa.bootstrpOneSampleTtest(diff(data), ...
    10000, 'MCPsol', 'none')

data = pse_table_prepost_eeg{2:3,4:end};
[p_uncorr,h_uncorr,~,~,obsStat] = mvpa.bootstrpOneSampleTtest(diff(data), ...
    10000, 'MCPsol', 'none')