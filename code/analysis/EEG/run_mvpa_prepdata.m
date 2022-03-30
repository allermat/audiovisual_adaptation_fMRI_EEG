% Input parameters 
trainMethod = 'sample-wise-avg';

subID = {'149','336','340','345','346'};

for i = 1:size(subID,2)
    analysisDir = fullfile(get_path('project'),'data',subID{i},'EEG','MVPA',trainMethod);
    if ~exist(analysisDir,'dir')
        mkdir(analysisDir);
    end
    
    I = struct();
    I.dir_analysis = analysisDir;
    % I.dir_preproc = fullfile(get_path('project'),'data',subID{i},'EEG','preproc_data','MVPA');
    I.dir_preproc = fullfile(get_path('project'),'data',subID{i},'EEG','MVPA','sample-wise');
    I.subID = subID{i};
    I.tr_method = trainMethod;
    mvpa.prepdata(I);
end
