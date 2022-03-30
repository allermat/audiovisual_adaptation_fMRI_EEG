%% Checking setup
[~,setupID] = system('hostname');
setupID = regexp(setupID,'[a-zA-Z_0-9-]*','match');
setupID = setupID{1};

if ~isempty(regexp(setupID,'^bb.+','once'))
    isServer = true;
else
    isServer = false;
    dbstop if error;
end

%% Opening parallel pool. 
% if there is no parallel pool running, open one. 
currPool = gcp('nocreate');
if isempty(currPool)
    if isServer
        parpool('local',16);
    else
        parpool('local');
    end
end

%% Input parameters 
subID = {'149','336','340','345','346'};
[feat,info] = deal(cell(size(subID)));
trainMethod = 'sample-wise-avg';

% Whether to display progress monitor
if isServer, progrMonitor = false; else progrMonitor = true; end

for iSubj = 1:numel(subID)
    
    dataFileName = [subID{iSubj},'_sw-avg_data.mat'];
    
    I = struct();
    I.dir_analysis = fullfile(get_path('project'),'data', ...
                              subID{iSubj},'EEG','PCM',trainMethod);
    if ~exist(I.dir_analysis,'dir')
        mkdir(I.dir_analysis);
    end
    I.dir_dataFile = fullfile(get_path('project'),'data',subID{iSubj},...
                              'EEG','PCM',trainMethod,dataFileName);
    I.condSel = 'nCatch';
    I.rdmLayout = 'pre-postL-postR';
    I.saveData = true;
    I.tr_method = trainMethod;
    I.timePoints = num2cell(-100:5:500);
    
    [feat{iSubj},info{iSubj}] = pcm_extractdata(I);
    
end

if ~isServer, dbclear if error; end