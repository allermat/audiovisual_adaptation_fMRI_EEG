function varargout = svmTrain(I)
% Method for training a support vector machine model on timeseries data
%
% USAGE:
%   trainS = svmTrain(I)
%   [trainS,filePath] = svmTrain(I)
% INPUT: 
%   I (structure): input settings with the following fields:  
%       cv_scheme: crossvalidation scheme, valid schemes: 'kf','loso','loxo'
%       dir_analysis: full path of the analysis directory. 
%       dir_dataFile: full path of the data file. 
%       dir_utils: full path of the utilities directory. 
%       k: number of folds in the k-fold crossvalidation 
%       nKFrep: number of repetitions of the k-fold crossvalidation
%       params: 1x2 cell array containing the vectors of possible parameters to
%           check during the gridsearch. 
%       subID:
%       svm_type: cc - C-SVC, nc - nu-SVC, er - epsilon-SVR, nr - nu-SVR
%       tr_cond: hyphen separated list of training conditions
%       tr_label: the label for training, valid labels: 'hand','resp',
%           'cond','sAhatInd','sVhatInd','sHatComm','sAhatANDsVhat'
%       tr_method: training method, valid methods: 
%           'sample-wise','sample-wise-beta','sample-wise-beta-rl',...
%           'sample-wise-bp','sample-wise-bp-avg','sample-wise-rl',...
%           'sample-wise-sc','sample-wise-mult','sample-wise-sm',...
%           'sample-wise-tf'
%       tr_timePoints: cell array of samples in ms to be used. 
%           If there are multiple time points specified in one cell, than the
%           features belonging to these time points are going to be treated as 
%           one set. 
%       progrMonitor (logical): wether to display a progress monitor 
% OUTPUT:  
%   filePath: full path to the output mvpa result structure saved on disc 
%       as a .mat file. The file contains the variables/fields below. 
%       The dimensions of the variables varies depending on the 
%       cross-validation method. 
%       
%       K-fold CV:
%       ==========
%       info: 1x1 struct with various information
%           
%       tr_examples: table of information about the examples used for 
%           training. Number of rows = number of examples,  
%       tr_grids: 1x2 struct array with the following fields
%           pass, grids, param1, param2 (for SVRs), bestPar
%       tr_groupings: M x N x O array with 
%           M = number of examples, 
%           N = I.nKFrep
%           O = number of time points 
%       tr_models: M x N x O x P array with M,N determined by the maximum 
%           size of the matrices generated from the output models of 
%           svmtrain, O = number of models per time point, P = number of
%           time points.
%       tr_scParam: M x N x O x P array with M = number of features, 
%           N = 2 (mean and std for Z-scoring), O = number of models per 
%           sample, P = number of samples. 
%       
%       Leave One Session Out (LOSO) CV:
%       ================================
%       info: 1x1 struct with various information about training
%           
%       tr_examples: table of information about the examples used for 
%           training. Number of rows = number of examples,  
%       tr_grids: 1x2 struct array with the following fields: 
%           pass, grids, param1, param2 (for SVRs), bestPar
%       tr_groupings: M x N array with 
%           M = number of examples 
%           N = number of time points 
%       tr_models: M x N x O x P array with M,N determined by the maximum 
%           size of the matrices generated from the output models of 
%           svmtrain, O = number of models per time point, P = number of
%           time points.
%       tr_scParam: M x N x O x P array with M = number of features, 
%           N = 2 (mean and std for Z-scoring), O = number of models per 
%           time point, P = number of time points.
%
%   trainS: The above results contained in a structure. 
%   filePath: The path to the results file if saved on disk
%       (saveFile == true)

% Copyright(C) 2017, Mate Aller
% allermat@gmail.com

%% Parsing input
p = inputParser;

% Valid values for certain inputs. 
validBpFreq = {'d','t','a','b','gl','gh'};
validCVschemes = {'kf','loro','loxo'};
validSCmethods = {'one-f','one-e','Z-f','Z-e','Z-ef','Z-fe'};
validSVMtypes = {'cc','er','nc','nr'};
validTfTypes = {'p','ph','pph'};
validTrLabels = {'hand','aloc','vloc'};
validTrMethods = {'sample-wise','sample-wise-avg','sample-wise-bp',...
                  'sample-wise-bp-avg','sample-wise-tf','sample-wise-tf-avg'};
validBinRegFun = {'mean','min','max'};

% Custom validation functions for certain inputs. 
    function checkTrTimePoints(x)
        if ~iscell(x) || ~isrow(x)
            error('mvpa:svmTrain:invalidInput','tr_timePoints must be a cell array of one row!');
        elseif any(~cellfun(@isrow,x)) || numel(unique(cellfun(@numel,x))) > 1
            error('mvpa:svmTrain:invalidInput',...
                'All elements of tr_timePoints must be row vectors of equal length!');
        end
    end

    function checkTfFreq(x)
        if ~iscell(x) || ~isrow(x) || numel(x) > 2
            error('mvpa:svmTrain:invalidInput',...
                'tf_freq must be a cell array of one row, maximum two elements!');
        elseif any(~cellfun(@isrow,x))
            error('mvpa:svmTrain:invalidInput',...
                'All elements of tf_freq must be row vectors!');
        end
    end

% Defining inputs.

addParameter(p,'cv_scheme','',@(x)any(validatestring(x,validCVschemes)));
addParameter(p,'dataS',[],@(x)validateattributes(x,{'struct'},{'nonempty'}));
addParameter(p,'dir_analysis','',@(x)exist(x,'dir'));
addParameter(p,'dir_dataFile','',@(x)exist(x,'file'));
addParameter(p,'k',[],@(x)validateattributes(x,{'numeric'},...
    {'scalar','integer','positive'}));
addParameter(p,'nKFrep',[],@(x)validateattributes(x,{'numeric'},...
    {'scalar','integer','positive'}));
addParameter(p,'nRepsAvg',1,@(x)validateattributes(x,{'numeric'},...
    {'scalar','integer','positive'}));
addParameter(p,'nTrialsToAvg',[],@(x)validateattributes(x,{'numeric'},...
    {'scalar','integer','positive'}));
addParameter(p,'params',{},@(x)validateattributes(x,{'cell'},{'nrows',1}));
addParameter(p,'permuteLabels',false,@islogical);
addParameter(p,'progrMonitor',true,@islogical);
addParameter(p,'saveFile',true,@islogical);
addParameter(p,'sc_method','',@(x)any(validatestring(x,validSCmethods)));
addParameter(p,'subID','',@ischar);
addParameter(p,'svm_type','',@(x)any(validatestring(x,validSVMtypes)));
addParameter(p,'bp_freq','',@(x)any(validatestring(x,validBpFreq)));
addParameter(p,'tf_freq',{},@checkTfFreq);
addParameter(p,'tf_type','',@(x)any(validatestring(x,validTfTypes)));
addParameter(p,'tr_cond','',@(x)validateattributes(x,{'char'},{'nrows',1}));
addParameter(p,'tr_label','',@(x)any(validatestring(x,validTrLabels)));
addParameter(p,'tr_method','',@(x)any(validatestring(x,validTrMethods)));
addParameter(p,'tr_timePoints',[],@checkTrTimePoints);
addParameter(p,'timeBinRegFun','mean',@(x)any(validatestring(x,validBinRegFun)));
% Parsing inputs.
parse(p,I);

% Assigning inputs to variables. 
cvScheme = p.Results.cv_scheme;
dataS = p.Results.dataS;
dirAnalysis = p.Results.dir_analysis;
dirDataFile = p.Results.dir_dataFile;
k = p.Results.k;
nKFrep = p.Results.nKFrep;
nRepsAvg = p.Results.nRepsAvg;
nTrialsToAvg = p.Results.nTrialsToAvg;
params = p.Results.params;
permuteLabels = p.Results.permuteLabels;
progrMonitor = p.Results.progrMonitor;
saveFile = p.Results.saveFile;
scMethod = p.Results.sc_method;
subID = p.Results.subID;
svmType = p.Results.svm_type;
bpFreq = p.Results.bp_freq;
tfFreq = p.Results.tf_freq;
tfType = p.Results.tf_type;
trCond = p.Results.tr_cond;
trLabel = p.Results.tr_label;
trMethod = p.Results.tr_method;
trTimePoints = p.Results.tr_timePoints;
timeBinRegFun = p.Results.timeBinRegFun;

% Custom checkings after parsing
requiredVars = {'cv_scheme','dir_analysis','dataFile','sc_method',...
    'subID','svm_type','tr_cond','tr_label','tr_method','tr_timePoints'};
if any(ismember(requiredVars,p.UsingDefaults))
    error('mvpa:svmTrain:missingInput',...
        'All required parameters must be specified!');
end

requiredVarsKF = {'k','nKFrep'};
if strcmp(cvScheme,'kf') && any(ismember(requiredVarsKF,p.UsingDefaults))
    error('mvpa:svmTrain:missingInput','''k'' and ''nKFrep'' must be specified!');
end

if ismember(trMethod,{'sample-wise-tf','sample-wise-tf-avg'})
    if ismember('tf_type',p.UsingDefaults)
        error('mvpa:svmTrain:missingInput','tr_type must be specified!');
    end
end

requiredVarsBP = {'bp_freq'};
if strcmp(trMethod,'sample-wise-bp') && any(ismember(requiredVarsBP,p.UsingDefaults))
    error('mvpa:svmTrain:missingInput','''bp_freq'' must be specified!');
end

requiredVarsBP = {'nTrialsToAvg'};
if ~isempty(regexp(trMethod,'-avg$','once')) && any(ismember(requiredVarsBP,p.UsingDefaults))
    error('mvpa:svmTrain:missingInput','''nTrialsToAvg'' must be specified!');
end

[p,I] = deal([]); %#ok

%% Checking...
% matlab version
if verLessThan('matlab', '8.3.0.532')
    error('mvpa:svmTrain:notSupportedMatlabVersion',...
        'MATLAB 2014a or newer is required! ');
end

% parallel pool 
currPool = gcp('nocreate');
if isempty(currPool)
    error('mvpa:svmTrain:missingParallelPool',...
        'Please open a parallel pool in order to run the function in parallel mode!');
else
    poolSize = currPool.NumWorkers;
end

% libsvm 
if isempty(regexp(path,'libsvm','once'))
    error('mvpa:svmTrain:missingToolbox','libsvm not found in path!');
end 

% parfor_progress
if progrMonitor
    if isempty(regexp(path,'parfor_progress','once'))
        error('mvpa:svmTrain:missingToolbox','parfor_progress not found in path!');
    end
end

%% Setting hyperparameters. 
switch svmType
    case 'cc'
        % Number of parameters required for the SVM
        nParam = 1;
        % For C-SVC, C is the cost parameter controlling the errors,
        % Here, log2(C) is given instead of C for computational convenience.
        % log2(C)
        if isempty(params), params{1} = 0; end
    case 'er'
        % Number of parameters required for the SVM
        nParam = 2;
        % For epsilon-SVR, C is the cost parameter controlling the errors,
        % e controls the width of the 'tube' which is fitted to the data.
        % Here, log2(C) and log2(e) is given for computational convenience.
        if isempty(params)
            params{1} = 0; % log2(C) or C = 1
            params{2} = -3.3219; % log2(e) or e = 0.1
        elseif numel(params) < 2
            error('mvpa:svmTrain:invalidInput',...
                'Either both hpyerparameters should be specified or none.');
        end
    case 'nc'
        % Number of parameters required for the SVM
        nParam = 1;
        % For nu-SVC, nu is the the upper and lower limit of errors and number
        % of support vectors respectively.
        % nu
        if isempty(params), params{1} = 0.5; end
    case 'nr'
        % Number of parameters required for the SVM
        nParam = 2;
        % For nuSVR, C is the cost parameter controlling the errors,
        % nu is the the upper and lower limit of errors and number of
        % support vectors respectively. Here, log2(C) is given instead of C 
        % for computational convenience.
        if isempty(params)
            params{1} = 0; % log2(C) or C = 1
            params{2} = 0.5; % nu
        elseif numel(params) < 2
            error('mvpa:svmTrain:invalidInput',...
                'Either both hpyerparameters should be specified or none! ');
        end
end

% Determining if gridsiearch is necessary
if nParam == 1
    if isscalar(params{1})
        doGridSearch = false;
    else
        doGridSearch = true;
        % Number of steps for the fine pass grid search.
        nStepsFine = 10;
    end
elseif nParam == 2
    if isscalar(params{1}) && isscalar(params{2})
        doGridSearch = false;
    else
        doGridSearch = true;
        % Number of steps for the fine pass grid search.
        nStepsFine = 10;
    end
end

% Starting the timer
cStartFun = clock;
ds = datestr(now);
fprintf('\nsvmtrain%s%s\n%s\n',repmat(' ',1,72-length(ds)),ds,repmat('-',1,80));

% Loading data and definig important variables
if isempty(dataS)
    dataS = load(dirDataFile);
end
% dataM = matfile(dirDataFile);
condDef = dataS.condDef;
misc = dataS.misc;

% Default value for trFreq if applicable
if ismember(trMethod,{'sample-wise-tf','sample-wise-tf-avg'})
    if isempty(tfFreq)
        if ismember(tfType,{'p','ph'})
            tfFreq = {misc.frequencies};
        else
            tfFreq = repmat({misc.frequencies},1,2);
        end
    end
end

% Converting training timepoints from milliseconds to seconds
trTimePointsSec = cellfun(@rdivide,trTimePoints,...
                              repmat({1000},size(trTimePoints)),'UniformOutput',false);

% Pre allocating arrays for collecting data
[avgGrouping,cvGroupings,models,acc,scParam] = deal(cell(nRepsAvg,1));
for iReps = 1:nRepsAvg
    
    % Choosing the examples for training accroding to the defined conditions
    % trExamples_feat is a 3D matrix with 
    % size(1) = nFeatures
    % size(2) = nTimePoints
    % size(3) = nExamples
    % Feature scaling is also done here within selectdata. 
    if ismember(trMethod,{'sample-wise-tf','sample-wise-tf-avg'})

        [trExamples_info,trExamples_feat,trIsExample,avgGrouping{iReps}] = ...
            mvpa.selectdata(dataS,trCond,trMethod,trLabel,trTimePointsSec,...
                            'tfFreq',tfFreq,'tfType',tfType,'progrMonitor',progrMonitor);
        
    else
        
        [trExamples_info,trExamples_feat,trIsExample,avgGrouping{iReps}] = ...
            mvpa.selectdata(dataS,trCond,trMethod,trLabel,trTimePointsSec,...
                            'progrMonitor',progrMonitor);
        
    end
    
    [nFeatures,nTimePoints,nExamples] = size(trExamples_feat);
    
    
    %% Preparing variables for cross-validation
    if strcmp(cvScheme,'kf')
        % K-fold cross-validation

        % Preparing groupings for the cross-validation
        misc.k = k;
        misc.nKFrep = nKFrep;
        [actCVgroupings,sRand] = mvpa.assigngroupings(cvScheme,trExamples_info,misc);
        % Cross-validation type specific input for the cross-validation 
        % function
        cvDetails.nKFrep = nKFrep;
        cvDetails.k = k;
        nModelsPerTimePoint = size(unique(actCVgroupings(:,1)),1)*size(actCVgroupings,2);
        % The maximum possible number of support vectors over the folds is the
        % number of examples in the largest of the training sets.
        maxnSVs = nExamples - minnival(actCVgroupings(:,1));
        % Choosing function according to the cross-validation scheme
        funCV = @mvpa.dokfoldcv;
        
    elseif strcmp(cvScheme,'loro')
        % Leave one run out cross-validation
        % Preparing groupings for the cross-validation
        actCVgroupings = mvpa.assigngroupings(cvScheme,trExamples_info,misc);
        % Cross-validation type specific input for the cross-validation 
        % function
        nModelsPerTimePoint = size(unique(actCVgroupings(:,1)),1);
        % The maximum possible number of support vectors over the folds is the
        % number of examples in the largest of the training sets.
        maxnSVs = nExamples - minnival(actCVgroupings(:,1));
        % Choosing function according to the cross-validation scheme
        funCV = @mvpa.dolosocv;
        
    elseif strcmp(cvScheme,'loxo')
        % Leave one example out cross-validation
        error('mvpa:svmTrain:missingFunctionality',...
              'Cross-validation scheme ''loxo'' is not yet implemented! ');
    end

    cvGroupings{iReps} = actCVgroupings;
    
    % General input for the crossvalidation function
    cvDetails.doGridSearch = doGridSearch;
    cvDetails.nFeatures = nFeatures;
    cvDetails.maxnSVs = maxnSVs;
    cvDetails.nModelsPerTimePoint = nModelsPerTimePoint;
    cvDetails.nParam = nParam;
    cvDetails.params = params;
    cvDetails.scMethod = scMethod;
    cvDetails.svmType = svmType;
    
    % Initializing arrays for data collection (padded with NaNs).
    modelsActRep = NaN(nFeatures+3,maxnSVs,nModelsPerTimePoint,nTimePoints);
    scParamActRep = NaN(nFeatures,2,nModelsPerTimePoint,nTimePoints);
    accActRep = NaN(3,nModelsPerTimePoint,nTimePoints);
    if doGridSearch
        [gridsCoarse,gridsFine,paramFine1,paramFine2,bestParCoarse,bestParFine] = deal(cell(nTimePoints,1));
        cvDetails.nStepsFine = nStepsFine;
    end

    %% Training models
    % Initializing progress monitor
    cStartTr = clock;
    fprintf('Training... \n');
    if progrMonitor
        parfor_progress(nTimePoints);
    end

    % Assigning training labels
    if ~ismember(trLabel,trExamples_info.Properties.VariableNames)
        error('mvpa:svmTrain:invalidVariableValue',...
              'The specified training label is not part of the dataset!');
    else
        trLabs = trExamples_info.(trLabel);
    end
    % Checking training labels
    if any(isnan(trLabs))
        error('mvpa:svmTrain:invalidVariableValue',...
              'Training labels can''t be NaNs!');
    elseif numel(unique(trLabs)) > 2 && any(ismember(svmType,{'cc','nc'}))
        error('mvpa:svmTrain:invalidVariableValue',...
              'The training label must have two levels for classification!');
    end
    
    % If permuted labels are required
    if permuteLabels
        if size(actCVgroupings,2) > 1
            error('mpva:svmTrain:missingFunctionality',...
                  ['Label permutation is not working ' ...
                   'if nKfReps > 1']);
        end
        rng('shuffle');
        for iFold = unique(actCVgroupings)'
            sel = find(actCVgroupings == iFold);
            trLabs(sel) = trLabs(sel(randperm(size(sel,1))));
        end
    end
    
    % If the number of time pointsis greater or equal to the parallel pool 
    % size, the parfor loop runs across the timepoints. 
    if nTimePoints >= poolSize
        
        parfor iTimePoints = 1:nTimePoints
            % Extracting the actual time sample. This yields a matrix with number
            % of rows = number of examples, number of columns = number of features.
            actTimePointFeats = squeeze(trExamples_feat(:,iTimePoints,:))';
            actTimePointLabs = trLabs;
            % Collecting data.
            [modelsActRep(:,:,:,iTimePoints),scParamActRep(:,:,:,iTimePoints),...
             accActRep(:,:,iTimePoints),cvMisc] = ...
                funCV(actTimePointFeats,actTimePointLabs,actCVgroupings,cvDetails);
            
            if doGridSearch
                gridsCoarse{iTimePoints} = cvMisc.gridsCoarse;
                gridsFine{iTimePoints} = cvMisc.gridsFine;
                paramFine1{iTimePoints} = cvMisc.paramFine1;
                if nParam == 2
                    paramFine2{iTimePoints} = cvMisc.paramFine2;
                end
                bestParCoarse{iTimePoints} = cvMisc.bestParCoarse;
                bestParFine{iTimePoints} = cvMisc.bestParFine;
            end
            % Advancing Progress monitor if applicable
            if progrMonitor
                parfor_progress;
            end
        end
        
        
        % If the number of time points is smaller than the parallel pool size, the 
        % parfor loop runs across the parameter space of the gridsearch. 
    else
        
        for iTimePoints = 1:nTimePoints
            
            % Extracting the actual time sample. This yields a matrix with number
            % of rows = number of examples, number of columns = number of features.
            actTimePointFeats = squeeze(trExamples_feat(:,iTimePoints,:))';
            actTimePointLabs = trLabs;
            % Collecting data.
            [modelsActRep{iReps}(:,:,:,iTimePoints),scParamActRep{iReps}(:,:,:,iTimePoints),...
             accActRep{iReps}(:,:,iTimePoints),cvMisc] = ...
                funCV(actTimePointFeats,actTimePointLabs,actCVgroupings,cvDetails);
            
            if doGridSearch
                gridsCoarse{iTimePoints} = cvMisc.gridsCoarse;
                gridsFine{iTimePoints} = cvMisc.gridsFine;
                paramFine1{iTimePoints} = cvMisc.paramFine1;
                if nParam == 2
                    paramFine2{iTimePoints} = cvMisc.paramFine2;
                end
                bestParCoarse{iTimePoints} = cvMisc.bestParCoarse;
                bestParFine{iTimePoints} = cvMisc.bestParFine;
            end
            % Advancing Progress monitor if applicable
            if progrMonitor
                parfor_progress;
            end
            
        end
        
    end
    
    models{iReps} = modelsActRep;
    scParam{iReps} = scParamActRep;
    acc{iReps} = accActRep;
    
    % Finalizing progress monitor.
    if progrMonitor
        parfor_progress(0);
    end
    fprintf('Training elapsed time (days hours:minutes:seconds) %s \n\n',...
            datestr(etime(clock,cStartTr)/86400,'dd HH:MM:SS'));

end

%% Organizing and saving the output 
trainS = struct();
trainS.info.cv_scheme = cvScheme;
if strcmp(cvScheme,'kf')
    trainS.info.k = k;
    trainS.info.nKFrep = nKFrep;
end
trainS.info.params = params;
trainS.info.sc_method = scMethod;
trainS.info.subID = subID;
trainS.info.svm_type = svmType;
if ismember(trMethod,{'sample-wise-tf','sample-wise-tf-avg'})
    trainS.info.tf_freq = tfFreq;
    trainS.info.tf_type = tfType;
elseif ismember(trMethod,{'sample-wise-bp','sample-wise-bp-avg'})
    trainS.info.bp_freq = bpFreq;
end
trainS.info.tr_cond = trCond;
trainS.info.tr_method = trMethod;
trainS.info.tr_timePoints = trTimePointsSec;
if any(cellfun(@numel,trTimePointsSec) > 1)
    trainS.info.timeBinRegFun = timeBinRegFun;
end
trainS.info.tr_elapsedTime = etime(clock,cStartFun);
trainS.info.tr_isExample = trIsExample;
if strcmp(cvScheme,'kf')
    trainS.info.tr_sRand = sRand;
end
[~,name,ext] = fileparts(dirDataFile);
trainS.info.tr_data_file = [name,ext];
trainS.info.tr_data_fs = misc.fs;
trainS.info.tr_label = trLabel;
if regexp(trMethod,'-avg$','once')
    trainS.info.tr_nTrialsToAvg = nTrialsToAvg;
    trainS.info.tr_avgGrouping = avgGrouping;
    trainS.info.tr_nRepsAvg = nRepsAvg;
end
trainS.info.sc_method = scMethod;

trainS.tr_examples = trExamples_info;

if doGridSearch && nRepsAvg == 1
    trainS.tr_grids(1).pass = 'coarse';
    trainS.tr_grids(1).grids = cat(ndims(gridsCoarse{1})+1,gridsCoarse{:});
    trainS.tr_grids(1).param1 = params{1}';
    if nParam == 2
        trainS.tr_grids(1).param2 = params{2}';
    end
    trainS.tr_grids(1).bestPar = cat(ndims(bestParCoarse{1})+1,bestParCoarse{:});
    trainS.tr_grids(2).pass = 'fine';
    trainS.tr_grids(2).grids = cat(ndims(gridsFine{1})+1,gridsFine{:});
    trainS.tr_grids(2).param1 = cat(ndims(paramFine1{1})+1,paramFine1{:});
    if nParam == 2
        trainS.tr_grids(2).param2 = cat(ndims(paramFine2{1})+1,paramFine2{:});
    end
    trainS.tr_grids(2).bestPar = cat(ndims(bestParFine{1})+1,bestParFine{:});
end

trainS.tr_groupings = cat(3,cvGroupings{:});

trainS.tr_models = cat(5,models{:});
models = [];
% Average over accuracies if multiple models were computed
trainS.tr_accuracies = nanmean(cat(4,acc{:}),4);

if ~strcmp(scMethod,'Z-f')
    trainS.tr_scParam = cat(5,scParam{:});
end

% Ordering the fields of the struct
trainS.info = orderfields(trainS.info);
trainS = orderfields(trainS); %#ok<*NASGU>

% Generating output file name.    
% Substring for the time samples
if size(trTimePoints,2) == 1
    if size(trTimePoints{1},2) > 1
        m = numel(trTimePoints{1});
        d = mean(diff(trTimePoints{1}));
        samplStr = sprintf('%d-m-%d-%d',mean(trTimePoints{1}),m,d);
    else
        samplStr = sprintf('%d',trTimePoints{1});
    end
elseif size(trTimePoints,2) > 1
    if any(cellfun(@numel,trTimePoints) > 1)
        m = numel(trTimePoints{1});
        d = mean(diff(trTimePoints{1}));
        trTimePoints = cellfun(str2func(timeBinRegFun),trTimePoints);
        samplStr = sprintf('%d-%d-%d-m-%d-%d',min(trTimePoints),...
            round(mean(diff(trTimePoints))),max(trTimePoints),m,d);
    else
        trTimePoints = cell2mat(trTimePoints);
        samplStr = sprintf('%d-%d-%d',min(trTimePoints),round(mean(diff(trTimePoints))),...
            max(trTimePoints));
    end
end
% Substring for training label
trLabelStr = trLabel;
 
% Substring for the frequencies
if ismember(trMethod,{'sample-wise-tf','sample-wise-tf-avg'})
    
   if ismember(tfType,{'p','ph'})
       % s = cell(size(tfFreq{1}));
       % for i = 1:size(tfFreq{1})
           % if mod(tfFreq{1}(i),1), s{i} = '-%.1f'; else s{i} = '-%d'; end
       % end
       % freqStr = strrep(sprintf(['_%s' [s{:}]],tfType,round(tfFreq{1})),'.',',');
       freqStr = '';
   else
       s = {cell(size(tfFreq{1})),cell(size(tfFreq{2}))};
       for i = 1:size(s,2)
           for j = 1:size(tfFreq{i},2)
               if mod(tfFreq{i}(j),1), s{i}{j} = '-%.1f'; else s{i}{j} = '-%d'; end
           end
       end
       freqStr = strrep(sprintf(['_f',[s{1}{:}],'_ph',[s{2}{:}]],tfFreq{1},tfFreq{2}),'.',',');      
   end

else
    freqStr = '';
end

fileName = sprintf('%s_tr-%s_%s_%s%s_%s.mat',svmType,trCond,trLabelStr,...
                   samplStr,freqStr,datestr(now,'yymmddHHMMSS'));
filePath = fullfile(dirAnalysis,fileName);
trainS.info.tr_filePath = filePath;
varargout{1} = trainS;

if saveFile
    % Saving outputs. 
    save(filePath,'-struct','trainS','-v7.3');
    varargout{2} = filePath;
end

end
