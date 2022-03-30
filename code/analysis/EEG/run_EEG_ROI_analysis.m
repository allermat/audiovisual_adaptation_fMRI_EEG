% Loading mvpa dataset
trainMethod = 'sample-wise-avg';
fileName = {'gr_tr-pre_gen-post.mat','gr_tr-pre_gen-pre.mat'};

for i = 1:numel(fileName)
    filePath = fullfile(get_path('project'),'data','group','EEG','MVPA',trainMethod,fileName{i});
    m = mvpares(filePath);
    
    % Determining time window for analysis
    time = m.getTrTimePoints;
    timeIdx = time >= 0.07 & time <= 0.13;
    
    % Loading individual objects
    srcFiles = m.getInfo.sourceFiles;
    srcFiles = cellfun(@repairPath,srcFiles,'UniformOutput',false);
    objList = cellfun(@mvpares,srcFiles,'UniformOutput',false);
    % Selecting data and performing statistical test
    if strcmp(fileName{i},'gr_tr-pre_gen-post.mat')
        indivData = cellfun(@getRecalIndex,objList,...
        repmat({'genTime'},size(objList)),...
        repmat({'tr'},size(objList)));
        temp = [indivData.RIbin];
    else
        indivData = cellfun(@getGenPerfEstimates,objList,...
        repmat({'genTime'},size(objList)),...
        repmat({'tr'},size(objList)));
        temp = [indivData.rf_avgFolds];
    end
    data = mean(temp(timeIdx,:));
    
    [pu(i),hu(i),~,~,st(i)] = mvpa.bootstrpOneSampleTtest(data,10000,'MCPsol','none');
    
    dataMean(i) = mean(data);
    dataSEM(i) = std(data)/sqrt(size(data,2));
end
eeg_roi_stats_table = table(pu',hu',st',dataMean',dataSEM');
eeg_roi_stats_table.Properties.VariableNames = {'p','h','t','mean','SEM'};
eeg_roi_stats_table.Properties.RowNames = {'RI','Acc'};
% Save results table
fileName = 'eeg_roi_stats.mat';
filePath = fullfile(get_path('project'),'data','group','EEG','MVPA',trainMethod,fileName);
save(filePath,'eeg_roi_stats_table');