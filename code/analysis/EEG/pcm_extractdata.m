function [feat_avg,info_avg] = pcm_extractdata(I)
% Function for extracting EEG data for PCM
%
% USAGE: 
%   rdmts = compRDMtimeseries(I)
% INPUT:
%   I (structure): input settings. Valid fields are: 
%       tr_method: 
%
% OUTPUT:
%   rdmts = RDM timeseries represented as a wrapped RDM set
%       (as implemented in the rsatoolbox)
%

% Copyright(C) 2017, Mate Aller
% allermat@gmail.com

% Parsing input
p = inputParser;

validTrMethods = {'sample-wise','sample-wise-avg'};
validBinRegFun = {'mean','min','max'};
validCVtypes = {'k-fold'};
validOutput = {'mean','byRun','glmByRun','glmPooled','rand'};

% Custom validation functions for certain inputs. 
    function checkTrTimePoints(x)
        if ~iscell(x) || ~isrow(x)
            error('mvpa:svmTrain:invalidInput','tr_timePoints must be a cell array of one row!');
        elseif any(~cellfun(@isrow,x)) || numel(unique(cellfun(@numel,x))) > 1
            error('mvpa:svmTrain:invalidInput',...
                'All elements of tr_timePoints must be row vectors of equal length!');
        end
    end

addParameter(p,'dir_analysis','',@(x)exist(x,'dir'));
addParameter(p,'dir_dataFile','',@(x)exist(x,'file'));
addParameter(p,'condSel','',@(x)validateattributes(x,{'char'},{'nrows',1}));
addParameter(p,'cvType','k-fold',@(x)any(validatestring(x,validCVtypes)));
addParameter(p,'k',4,@(x)validateattributes(x,{'numeric'},...
                                            {'scalar','positive','integer'}));
addParameter(p,'rdmLayout','pre',@(x)validateattributes(x, ...
                                                  {'char'},{'nrows',1}));
addParameter(p,'tr_method','',@(x)any(validatestring(x,validTrMethods)));
addParameter(p,'timePoints',[],@checkTrTimePoints);
addParameter(p,'timeBinRegFun','mean',@(x)any(validatestring(x,validBinRegFun)));
addParameter(p,'output','mean',@(x) ismember(x,validOutput));

% Parsing inputs.
parse(p,I);

% Assigning inputs to variables. 
dirAnalysis = p.Results.dir_analysis;
dirDataFile = p.Results.dir_dataFile;
condSel = p.Results.condSel;
cvType = p.Results.cvType; %#ok
k = p.Results.k;
rdmLayout = p.Results.rdmLayout;
trMethod = p.Results.tr_method;
timePoints = p.Results.timePoints;
timeBinRegFun = str2func(p.Results.timeBinRegFun);
output = p.Results.output;

% Loading data as prepared for mvpa
dataM = matfile(dirDataFile);
% Extracting relevant fields
info = dataM.info;
condDef = dataM.condDef;
% Converting time input to seconds
timePointsSec = cellfun(@rdivide,timePoints,...
                        repmat({1000},size(timePoints)),...
                        'UniformOutput',false);

% Extracting features
trialSel = mvpa.selectexamples(condSel,condDef,info);
info = info(trialSel,:);
info = ersa.createGroupingVariable(info,rdmLayout);
info.idx = (1:size(info,1))';

feat = ersa.extractfeatures(dataM,trialSel,trMethod,timePointsSec);

% Computing RDMs for each selected timepoint separately
nConds = size(unique(info.rdmGrVar),1);

if strcmp(output,'mean')
    trialIdxToAvg = varfun(@(x) x,info,'InputVariables',{'idx'}, ...
        'GroupingVariables',{'rdmGrVar'},...
        'OutputFormat','cell');
    [~,ia] = unique(info.rdmGrVar);
elseif strcmp(output,'rand')
    nTrialsToAvg = 20;
    nRuns = size(unique(info(:,{'session','run'})),1);
    info = grouptrials(info,condDef,nTrialsToAvg);
    [c,ia] = unique(info.avgExampleID);
    ia = ia(~isnan(c));
    trialIdxToAvg = varfun(@(x) x,info,'InputVariables',{'idx'}, ...
        'GroupingVariables',{'avgExampleID'},'OutputFormat','cell');
    info = removevars(info,'avgExampleID');
elseif ismember(output,{'byRun','glmByRun'})
    info.partVec = info.run;
    u = unique(info(:,{'session','run'}));
    for i = 1:size(u,1)
        info.partVec(ismember(info(:,{'session','run'}),u(i,:),'rows')) = i;
    end
    trialIdxToAvg = varfun(@(x) x,info,'InputVariables',{'idx'}, ...
        'GroupingVariables',{'partVec','rdmGrVar'},...
        'OutputFormat','cell');
    [~,ia] = unique(info(:,{'partVec','rdmGrVar'}),'rows');
    % Remove conditions which do not have at least 20 trials
    nTrialsMin = 20;
    toReject = cellfun(@numel,trialIdxToAvg) < nTrialsMin;
    trialIdxToAvg(toReject) = [];
    ia(toReject) = [];
    % Randomly select exactly 20 trials for each condition
    trialIdxToAvg = cellfun(@(x) x(randperm(nTrialsMin)),trialIdxToAvg,...
                          'UniformOutput',false);
end

if ismember(output,{'mean','byRun','rand'})
    % Averaging over examples for the same condition
    feat_avg = cellfun(@(x) mean(feat(:,:,x),3),trialIdxToAvg, ...
        'UniformOutput',false);
    feat_avg = cat(3,feat_avg{:});
    % Creating info table for average features, leaving session specific info
    % out
    info_avg = info(ia,6:end);
    if strcmp(output,'rand')
        temp = mod(randperm(size(info_avg,1))',nRuns);
        temp(temp == 0) = nRuns;
        info_avg.partVec = temp;
    end
elseif strcmp(output,'glmByRun')
    % First get rid of trials which are not needed
    tempInfo = info(sort(cat(1,trialIdxToAvg{:})),:);
    tempFeat = feat(:,:,sort(cat(1,trialIdxToAvg{:})),:);
    % Update trialIdxToAvg to reflect the selected data
    tempInfo.idx = (1:size(tempInfo,1))';
    tempTrlIdxToAvg = varfun(@(x) x,tempInfo,'InputVariables',{'idx'}, ...
        'GroupingVariables',{'partVec','rdmGrVar'},...
        'OutputFormat','cell');
    [~,tempIa] = unique(tempInfo(:,{'partVec','rdmGrVar'}),'rows');
    % Generate design matrix
    partitions = unique(tempInfo.partVec);
    % Overall mean
    X = [ones(size(tempInfo,1),1),zeros(size(tempInfo,1),numel(partitions))];
    % Run specific means
    for i = 1:numel(partitions)
        X(tempInfo.partVec == partitions(i),1+i) = 1;
    end
    % Condition specific means
    nColsX = size(X,2);
    X = cat(2,X,zeros(size(tempInfo,1),numel(tempTrlIdxToAvg)));
    for i = 1:numel(tempTrlIdxToAvg)
        X(tempTrlIdxToAvg{i},nColsX+i) = 1;
    end
    % Model responses
    s = size(tempFeat);
    betas = NaN(s(1),s(2),size(X,2));
    for iFeat = 1:s(1)
        for iTime = 1:s(2)
            y = squeeze(tempFeat(iFeat,iTime,:));
            betas(iFeat,iTime,:) = pinv(X)*y;
        end
    end
    feat_avg = betas(:,:,nColsX+1:end);
    info_avg = tempInfo(tempIa,6:end);
end

if ~isempty(dirAnalysis)
    % Generating output file name.
    % Substring for the time samples
    if size(timePoints,2) == 1
        if size(timePoints{1},2) > 1
            m = numel(timePoints{1});
            d = mean(diff(timePoints{1}));
            samplStr = sprintf('%d-m-%d-%d',mean(timePoints{1}),m,d);
        else
            samplStr = sprintf('%d',timePoints{1});
        end
    elseif size(timePoints,2) > 1
        if any(cellfun(@numel,timePoints) > 1)
            m = numel(timePoints{1});
            d = mean(diff(timePoints{1}));
            timePoints = cellfun(timeBinRegFun,timePoints);
            samplStr = sprintf('%d-%d-%d-m-%d-%d',min(timePoints),...
                round(mean(diff(timePoints))),max(timePoints),m,d);
        else
            timePoints = cell2mat(timePoints);
            samplStr = sprintf('%d-%d-%d',min(timePoints),round(mean(diff(timePoints))),...
                max(timePoints));
        end
    end

    fileName = sprintf('pcm_%s_%s_%s.mat',condSel,samplStr,datestr(now,'yymmddHHMMSS'));
    filePath = fullfile(dirAnalysis,fileName);
    
    % Saving outputs.
    save(filePath,'feat_avg','info_avg','I','-v7.3');
end
end

function [info,seed] = grouptrials(info,condDef,nTrialsToAvg)

blockTypes = unique(info.blocktype);
actConds = condDef.condition(ismember(condDef.blocktype,blockTypes));
info.avgExampleID = NaN(size(info,1),1);
nExamplesOverall = 0;
hands = unique(info.hand);
adaptDir = unique(info.adaptdir);
sessions = unique(info.session);
catchLevels = [true,false];
factors = fullfact([size(actConds,1),numel(adaptDir)]);
% factors = fullfact([size(actConds,1),numel(sessions),numel(catchLevels)]);

rng('shuffle');
seed = rng;

for iFact = 1:size(factors,1)
    % Making sure only trials with same response hand are
    % averaged 
    isExample = info.condition == actConds(factors(iFact,1)) & ...
                info.adaptdir == adaptDir(factors(iFact,2));
%               info.catch_trial == catchLevels(factors(iFact,3));
    if all(~isExample)
        warning('mvpa:averagetirals:grouptrials:missingExample',...
              'No example was found for condition %d, ', ...
              actConds(factors(iFact,1)));
    end

    % Number of examples after averaging in this condition
    nExamplesOut = floor(sum(isExample)/nTrialsToAvg);
    if nExamplesOut < 1
        warning('pcm_extractdata:grouptrials:missingExample',...
              'Not enough examples for averaging for condition %d, ', ...
              actConds(factors(iFact,1)));
    end
    % Removing the sruplus examples
    nSurplus = mod(sum(isExample),nTrialsToAvg);
    isExample(randsample(find(isExample),nSurplus)) = false;
    % Assigning example IDs to 
    avgExampleID = mod(randperm(sum(isExample))',nExamplesOut);
    avgExampleID(avgExampleID == 0) = nExamplesOut;
    % This makes sure that all averaged examples have a unique ID
    info.avgExampleID(isExample) = avgExampleID + nExamplesOverall;
    nExamplesOverall = nExamplesOverall + nExamplesOut;
end

end
