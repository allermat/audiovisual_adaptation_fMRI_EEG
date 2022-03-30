% Figure 3
clearvars;
close all;
% Panel A: fMRI decoding accuracy and recalibration index
dataPath = fullfile(get_path('project'),'results','data','mvpa_fmri.mat');
if ~exist(dataPath,'file')
    out = run_mvpa_fmri;
    save(dataPath,'out');
end
[fmri_mvpa_Acc_table, fmri_mvpa_RI_table] = figure_fmri_mvpa_index(dataPath);

% Calculate bootstrap-based one sample ttests
data = fmri_mvpa_Acc_table{:,4:end};
[p_uncorr(:,1),h_uncorr,~,~,obsStat(:,1)] = mvpa.bootstrpOneSampleTtest(...
    data, 10000, 'MCPsol', 'none')

data = fmri_mvpa_RI_table{:,4:end};
[p_uncorr(:,2),h_uncorr,~,~,obsStat(:,2)] = mvpa.bootstrpOneSampleTtest(...
    data, 10000, 'MCPsol', 'none')

% For multiple comparison correction across ROIS, we use Benjamini-Hochberg 
% algorithm which is taken from Matlab Central Fie Exchange
% (https://uk.mathworks.com/matlabcentral/fileexchange/27418-fdr_bh)
[h_corr,crit_p,adj_ci_cvrg,p_corr] = fdr_bh(p_uncorr(:,1), 0.05)

[h_corr,crit_p,adj_ci_cvrg,p_corr] = fdr_bh(p_uncorr(:,2), 0.05)

% Panel B: fMRI neurometric functions
[fmri_NF_table, fmri_NF_data_table] = figure_fmri_NF(dataPath);

% Panel C: fMRI representational dissimilarity matrices and
% multidimensional scaling
figures_MDS;