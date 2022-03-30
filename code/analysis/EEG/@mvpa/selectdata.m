function [info,feat,isExampleOrig,trialGrouping] = selectdata(dataM,cond,trMethod,trLabel,timePoints,varargin)
% Selects all necessary data (features, labels) for mvpa
% 
% INPUT:
%   dataM (matfile object): matfile object referencing to the data
%       for mvpa. It should contain the following fields: feat,
%       info, misc, condDef
%   cond (string): Condition name for data to be selected
%   trMethod (string): training method name
%   trLabel (string): training label name
%   timePoints (cell array): cell array of time points to be
%       selected
%   
% OUTPUT:
%   info (table): table of nExamples rows, containing information
%       corresponding the features in feat
%   feat (numerical array): features array of size nFeatures x
%       nTimePoints x nExamples
%   isExample (logical array): logical vector indicating the
%       selected examples in the original dataset. 
%   trialGrouping (numeric array): logical vector indicating the
%       grouping of the original examples for averaged examples
%       (only aplicable if trMethod is averaged).
%

% Parsing input
p = inputParser;

% Valid values for certain inputs.
validTrMethods = {'sample-wise','sample-wise-avg','sample-wise-tf',...
                  'sample-wise-tf-avg'};
validTrLabels = {'hand','aloc','vloc'};
validTfTypes = {'p','ph','pph'};
% Custom validation functions for certain inputs. 
    function checkTimePoints(x)
        if ~iscell(x) || ~isrow(x)
            error('mvpa:extractfeatures:invalidInput',...
                'tr_samples must be a cell array of one row!');
        elseif any(~cellfun(@isrow,x)) || numel(unique(cellfun(@numel,x))) > 1
            error('mvpa:extractfeatures:invalidInput',...
                'All elements of tr_timePoints must be row vectors of equal length!');
        end
    end

% Defining inputs.
addRequired(p,'dataM',@(x)validateattributes(x,{'struct','matlab.io.MatFile'},{'nonempty'}));
addRequired(p,'cond',@(x)validateattributes(x,{'char'},{'nrows',1}));
addRequired(p,'trMethod',@(x)any(validatestring(x,validTrMethods)));
addRequired(p,'trLabel',@(x)any(validatestring(x,validTrLabels)));
addRequired(p,'timePoints',@checkTimePoints);
addParameter(p,'progrMonitor',true,@islogical);
addParameter(p,'tfFreq',{},@iscell);
addParameter(p,'tfType','',@(x)any(validatestring(x,validTfTypes)));
addParameter(p,'trialGrouping',[],@(x)validateattributes(x,{'numeric'},{'column'}));

% Parsing inputs.
parse(p,dataM,cond,trMethod,trLabel,timePoints,varargin{:});

% Assigning inputs to variables.
dataM = p.Results.dataM;
cond = p.Results.cond;
trMethod = p.Results.trMethod;
trLabel = p.Results.trLabel;
timePoints = p.Results.timePoints;
progrMonitor = p.Results.progrMonitor;
tfFreq = p.Results.tfFreq;
tfType = p.Results.tfType;
trialGrouping = p.Results.trialGrouping;

% Custom checkings after parsing
if ismember(trMethod,{'sample-wise-tf','sample-wise-tf-avg'}) && ...
        any(cellfun(@isempty,{tfFreq,tfType}))
    error('mvpa:selectdata:missingInput',...
        'tfFreq and tfType must be specified!');
end

if strcmp(trMethod,'sample-wise-bp') && isempty(bpFreq)
    error('mvpa:selectdata:missingInput','bpFreq must be specified!');
end

% Selecting the info field name
if ismember(trMethod,{'sample-wise-beta','sample-wise-beta-rl','sample-wise-sm-beta'})
    if strcmp(trLabel,'stim')
        infoStr = 'info_stim';
    elseif strcmp(trLabel,'resp')
        infoStr = 'info_resp';
    elseif strcmp(trLabel,'hand')
        error('mvpa:selectdata:missingFunctionality',...
            'This analysis is not yet implemented!');
    end
else
    infoStr = 'info';
end
dataInfo = dataM.(infoStr);
condDef = dataM.condDef;
misc = dataM.misc;

% Selecting examples according to the input condition. 
isExampleOrig = mvpa.selectexamples(cond,condDef,dataInfo);
if sum(isExampleOrig) == 0
    error('mvpa:selectdata:requestedDataNotPresent',...
          'No examples of the specified condition were found in the dataset!');
end

if ismember(trMethod,{'sample-wise','sample-wise-avg'})
    
    featStr = 'feat';
    dataFeat = dataM.(featStr);
    if strcmp(trMethod,'sample-wise-avg')
        if isempty(trialGrouping)
            [tempInfo,dataInfo,dataFeat] = mvpa.averagetrials(dataInfo,dataFeat, ...
                                                         condDef,'rand',16);
            trialGrouping = tempInfo.avgExampleID;
            tempInfo = []; %#ok
        else
            [~,dataInfo,dataFeat] = mvpa.averagetrials(dataInfo,dataFeat, ...
                                                  condDef,'preset',16,...
                                                  trialGrouping);
        end
        % Because of the averaging the number of trials in the
        % info has changed, hence the examples need to be updated
        isExample = mvpa.selectexamples(cond,condDef,dataInfo);
    else
        isExample = isExampleOrig;
    end
    
    info = dataInfo(isExample,:);
    
    dataFeatSize = size(dataFeat);
    nFeatures = dataFeatSize(1)*size(timePoints{1},2);
    nSamples = size(timePoints,2);
    nExamples = sum(isExample);
    feat = NaN(nFeatures,nSamples,nExamples);
    
    % Initializing progress monitor
    cStartExtr = clock;
    fprintf('Selecting features...\n');
    if progrMonitor
        parfor_progress(size(timePoints,2));
    end
    
    for i = 1:size(timePoints,2)
        
        actSample = mvpa.indsamplecustom(timePoints{i},dataFeatSize(2),misc.fs,misc.timeOnset);
        % Error if any specified training sample is out of the boundaries of the
        % data.
        if any(isnan(actSample))
            error('mvpa:selectdata:invalidSample',...
                'At least one specified sample is not part of the data set!');
        end
        
        % Select the actual time sample(s) and the appropriate examples
        temp = dataFeat(:,actSample,isExample);
        
        % Reshaping the array to nFeatures x 1 x nExamples 
        s = size(temp);
        feat(:,i,:) = reshape(temp,s(1)*s(2),1,s(3));
        
        % Advancing progress bar
        if progrMonitor
            parfor_progress;
        end
        
    end
    
    
elseif ismember(trMethod,{'sample-wise-tf','sample-wise-tf-avg'})
    
    if ismember(tfType,{'p','ph'})
        
        if strcmp(tfType,'p')
            featStr = 'feat_p';
        else
            featStr = 'feat_ph';
        end
        fields = whos(dataM);
        fields = {fields.name}';
        if ~ismember(featStr,fields)
            error('mvpa:selectdata:missingFeatureType',...
                  'The specified feature is missing from the dataset.');
        end
        
        dataFeat = dataM.(featStr);
        if strcmp(trMethod,'sample-wise-avg')
            if isempty(trialGrouping)
                [tempInfo,dataInfo,dataFeat] = mvpa.averagetrials(dataInfo,dataFeat, ...
                                                             condDef,'rand',16);
                trialGrouping = tempInfo.avgExampleID;
                tempInfo = []; %#ok
            else
                [~,dataInfo,dataFeat] = mvpa.averagetrials(dataInfo,dataFeat, ...
                                                      condDef,'preset',16,...
                                                      trialGrouping);
            end
            % Because of the averaging the number of trials in the
            % info has changed, hence the examples need to be updated
            isExample = mvpa.selectexamples(cond,condDef,dataInfo);
        else
            isExample = isExampleOrig;
        end
        
        info = dataInfo(isExample,:);
        
        dataFeatSize = size(dataFeat);
        nFeatures = dataFeatSize(1)*size(timePoints{1},2)*size(tfFreq{1},2);
        nSamples = size(timePoints,2);
        nExamples = sum(isExample);
        feat = NaN(nFeatures,nSamples,nExamples);
        
        % Initializing progress monitor
        cStartExtr = clock;
        fprintf('Selecting features...\n');
        if progrMonitor
            parfor_progress(size(timePoints,2));
        end
        
        
        for i = 1:size(timePoints,2)
            
            actSample = mvpa.indsamplecustom(timePoints{i},dataFeatSize(3),misc.fs,misc.timeOnset);
            % Error if any specified training sample is out of the boundaries of the
            % data.
            if any(isnan(actSample))
                error('mvpa:selectdata:invalidSample',...
                    'At least one specified sample is not part of the data set!');
            end
            
            actFreq = tfFreq{1};
            % Error if any specified frequency is not part of the data
            if any(~ismember(actFreq,misc.frequencies))
                error('mvpa:selectdata:invalidFrequency',...
                    'At least one specified frequency is not part of the data set!');
            end
            actFreqInd = find(ismember(misc.frequencies,actFreq));
            temp = dataFeat(:,actFreqInd,actSample,isExample); %#ok<*FNDSB>
                        
            % 1. Reshaping the array to nFeatures x 1 x nExamples
            s = size(temp);
            feat(:,i,:) = reshape(temp,s(1)*s(2)*s(3),1,s(4));
            
            % Advancing progress bar
            if progrMonitor
                parfor_progress;
            end
            
        end

    end
    
end

% Finalizing progress monitor 
if progrMonitor
    parfor_progress(0);
end
fprintf('Feature extraction elapsed time (days hours:minutes:seconds) %s \n\n',...
    datestr(etime(clock,cStartExtr)/86400,'dd HH:MM:SS'));

end

