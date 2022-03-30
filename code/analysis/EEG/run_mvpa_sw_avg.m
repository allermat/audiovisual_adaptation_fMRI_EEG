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

trainMethod = 'sample-wise-avg';

% Settings for training
trSettings(1).cond = 'pre';
trSettings(1).label = 'aloc';
trSettings(1).fileName = '';
% Settings for generalizing
genSettings(1).cond = {'pre','nCatch'};
genSettings(1).genTime = {'tr','tr'};

% Whether to display progress monitor
if isServer, progrMonitor = false; else progrMonitor = false; end

for iSubj = 1:numel(subID)
    
    dataFileName = [subID{iSubj},'_sw-avg_data.mat'];
    
    for iTrain = 1:size(trSettings,2)
        
        if isempty(trSettings(iTrain).fileName)
            I = struct();
            I.cv_scheme = 'kf';
            I.dir_analysis = fullfile(get_path('project'),'data',subID{iSubj},'EEG','MVPA',trainMethod);
            I.dir_dataFile = fullfile(I.dir_analysis,dataFileName);
            I.k = 4;
            I.nKFrep = 1;
            I.nRepsAvg = 50;
            I.nTrialsToAvg = 16;
            % log2(C)
            % I.params{1} = -7:2;
            % log2(e)
            % I.params{2} = -5:3;
            I.progrMonitor = progrMonitor;
            I.saveFile = false;
            I.sc_method = 'Z-e';
            I.subID = subID{iSubj};
            I.svm_type = 'nr';
            I.tr_cond = trSettings(iTrain).cond;
            I.tr_label = trSettings(iTrain).label;
            I.tr_method = trainMethod;
            I.tr_timePoints = cellfun(@plus,num2cell(-50:5:500), ...
                              repmat({-50:5:0},1,size(-50:5:500,2)),'UniformOutput',false);
            I.timeBinRegFun = 'max';
            % I.tr_timePoints = num2cell(-100:5:600);
            
            trainS = mvpa.svmTrain(I);
        end
        
        if ~isempty(genSettings(iTrain))
            
            for iGen = 1:numel(genSettings(iTrain).cond)
                I = struct();
                if isempty(trSettings(iTrain).fileName)
                    % I.dir_trainFile = fn;
                    I.trainS = trainS;
                else
                    saveDf = cd(fullfile(get_path('project'),'data',subID{iSubj},'EEG','MVPA',trainMethod));
                    listing = dir;
                    fileNames = {listing.name}';
                    fileNames = regexp(fileNames,trSettings(iTrain).fileName,'match','once');
                    fileNames = fileNames(~strcmp(fileNames,''));
                    if size(fileNames,1) > 1
                        error('More than one file found with the specified match string!')
                    else
                        fileName = fileNames{1};
                    end
                    cd(saveDf);
                    I.dir_trainFile = fullfile(get_path('project'),'data',subID{iSubj},'EEG','MVPA',...
                                               trainMethod,fileName);
                end
                I.dir_dataFile = fullfile(get_path('project'),'data',subID{iSubj},'EEG','MVPA',...
                                          trainMethod,dataFileName);
                I.gen_cond = genSettings(iTrain).cond{iGen};
                I.gen_time = genSettings(iTrain).genTime{iGen};
                I.progrMonitor = progrMonitor;
                
                mvparesObj = mvpa.svmGeneralize(I);
                mvparesObj = mvpa.addGenPerfEstimates(mvparesObj);
                if strcmp({I.gen_cond},{'nCatch','all'})
                    mvparesObj = mvpa.addNeurometricFunctions(mvparesObj);
                    mvparesObj = mvpa.addRecalIndex(mvparesObj);
                end
                
            end
            
        end
        
    end
    
end

if ~isServer, dbclear if error; end