%% Loading grand average data
loc = [-12,-5,-2,0,2,5,12]';
pre_list = strcat('fteeg_gr_avg_ERP_pre-test_',cellfun(@num2str,num2cell(loc),'UniformOutput',false));
post_L_list = strcat('fteeg_gr_avg_ERP_post-test_L_',cellfun(@num2str,num2cell(loc),'UniformOutput',false));
post_R_list = strcat('fteeg_gr_avg_ERP_post-test_R_',cellfun(@num2str,num2cell(loc),'UniformOutput',false));

saveDir = cd(fullfile(get_path('project'),'results','data'));
pre_data = struct2cell(cellfun(@load,pre_list));
post_L_data = struct2cell(cellfun(@load,post_L_list));
post_R_data = struct2cell(cellfun(@load,post_R_list));
cd(saveDir);

%% Plotting figures
cfg = struct();
cfg.layout = 'acticap-64ch-standard2';
cfg.channel = {'FC1','FC2','Cz','CP1','CP2','C1','C2','CPz'};
cfg.xlim = [-0.05,0.5]; 
cfg.linewidth = 1.5;

% Pre, PostL and PostR ERPs for 0 location
condStr = {'pre 0','postL 0','postR 0'};
condData = {pre_data{loc == 0},post_L_data{loc == 0},post_R_data{loc == 0}};
figure(); ft_singleplotER(cfg,condData{:});
legend(condStr);
% save plot data to table
[temp_eeg_recal,temp_time_recal,temp_aloc_recal] = deal([]);
temp_cond_recal = {};
for i = 1:numel(condStr)
    actData = condData{i};
    chanIdx = find(ismember(actData.label,cfg.channel));
    timeIdx = find(actData.time >= cfg.xlim(1) & ...
                   actData.time <= (cfg.xlim(2)+eps)); % Adding a tiny bit to avoid floatin point precision issues
    temp_eeg_recal = cat(1,temp_eeg_recal,mean(actData.avg(chanIdx,timeIdx))');
    temp_time_recal = cat(1,temp_time_recal,actData.time(timeIdx)');
    ts = strsplit(condStr{i},' ');
    temp_cond_recal = cat(1,temp_cond_recal,...
                          repmat(ts(1),numel(timeIdx),1));
    temp_aloc_recal = cat(1,temp_aloc_recal,zeros(numel(timeIdx),1));
                          
end
eeg_erp_recal_table = table(temp_cond_recal,temp_aloc_recal,temp_time_recal,temp_eeg_recal,...
                        'VariableNames',{'cond','aloc','time','eeg'});

% Pretest ERPs across locations
cmap = viridis;
cfg.graphcolor = cmap(round(linspace(1,255,numel(loc))),:);
figure(); ft_singleplotER(cfg,pre_data{:});
legend(cellfun(@num2str,num2cell(loc),'UniformOutput',false));
% save plot data to table
[temp_eeg_pre,temp_time_pre,temp_aloc_pre] = deal([]);
temp_cond_pre = {};
for i = 1:numel(loc)
    actData = pre_data{i};
    chanIdx = find(ismember(actData.label,cfg.channel));
    timeIdx = find(actData.time >= cfg.xlim(1) & ...
                   actData.time <= (cfg.xlim(2)+eps));
    temp_eeg_pre = cat(1,temp_eeg_pre,mean(actData.avg(chanIdx,timeIdx))');
    temp_time_pre = cat(1,temp_time_pre,actData.time(timeIdx)');
    temp_cond_pre = cat(1,temp_cond_pre,repmat({'pre'},numel(timeIdx),1));
    temp_aloc_pre = cat(1,temp_aloc_pre,repmat(loc(i),numel(timeIdx),1));
                          
end
eeg_erp_pre_table = table(temp_cond_pre,temp_aloc_pre,temp_time_pre,temp_eeg_pre,...
                        'VariableNames',{'cond','aloc','time','eeg'});

% Loading mvpa dataset
trainMethod = 'sample-wise-avg';
fileName = 'gr_tr-pre_gen-post.mat';
filePath = fullfile(get_path('project'),'results','data',fileName);
m_post = mvpares(filePath);

fileName = 'gr_tr-pre_gen-pre.mat';
filePath = fullfile(get_path('project'),'results','data',fileName);
m_pre = mvpares(filePath);

% Loading ROI analysis stats
fileName = 'eeg_roi_stats.mat';
filePath = fullfile(get_path('project'),'results','data',fileName);
load(filePath,'eeg_roi_stats_table');

% Decoding accuracy and recal index
% This plots the recalibration index in percentages, the figure in the
% manuscript displays it in fractions
m_post.showRecalIndex('kind','bin');
m_pre.showGenPerf('perfIdx','rf');

% save plot data to table
temp = m_pre.getGenPerfEstimates;
temp_acc = temp.rf_avgFolds;
temp_acc_err = temp.rf_avgFolds_err;
temp_time = m_pre.getTrTimePoints;
eeg_decode_acc_table = table(temp_time,temp_acc,temp_acc_err,'VariableNames',...
                             {'time','mean','sem'});
temp = m_post.getRecalIndex;
temp_RI = temp.RIbin/100; % Converting to fractions for table
temp_RI_err = temp.RIbin_err/100; % Converting to fractions for table
temp_time = m_post.getTrTimePoints;
eeg_decode_RI_table = table(temp_time,temp_RI,temp_RI_err,'VariableNames',...
                            {'time','mean','sem'});