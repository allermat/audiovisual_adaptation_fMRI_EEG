function varargout = svmGeneralize(I)
% Method for generalizing a support vector machine model on timeseries data
% 
% USAGE:
%   [mvparesObj,filePath] = svmGeneralize(I)
%   [mvparesObj,genS] = svmGeneralize(I)
% INPUT: 
%   I (struct): input settings. 
%       Required fields:
%           dir_dataFile (string): path to the MVPA dataset
%           dir_trainFile (string): path to the trained mvpares dataset
%           gen_cond (string): generalization condition
%           gen_time (string): indicating the generalization time. Possible 
%               values: 'tr','tr_x_tr' for trainin time and training time 
%               by training time respectively
%       Optional fields:
%           customTag (string): arbitrary tag to be added to the name of
%               the output file
%           progrMonitor (logical): wether to display a progress monitor,
%              default: true
% OUTPUT:
%   The mvpa result structure contains the variables/fields of the
%       training mvpa result file and the following fields with the 
%       generalization results:
%       gen_examples: table of information about the examples used for 
%           generalization. Number of rows = number of examples, 
%       gen_accuracies: M x N x O x P array with 
%           M = 3 (accuracy,MSE,R2), 
%           N = number of models per sample
%           O = number of training samples
%           P = number of generalization samples per training sample
%       gen_predlabels: M x N x O x P array with 
%           M = maximum number of generalization examples during the 
%               cross-validation padded with NaNs if the actual number of 
%               test examples is smaller than that)
%           N = number of models per sample
%           O = number of samples 
%           P = number of generalization samples per training sample
%       The info variable is extended with further generalization
%           details.
% 
%   mvparesObj: The above results loaded as an mvpares object. 
%   filePath: path to the saved mvpares data structure 
%       (if saveFile == ture)
%   genS: mvpares data structure as a strucure array
%       (if saveFile == false)
% 

% Copyright(C) 2017, Mate Aller
% allermat@gmail.com

%% Parsing input, checking matlab
p = inputParser;

% Valid values for certain inputs.
validGenTimes = {'tr','tr_x_tr'};

% Defining inputs.
addParameter(p,'customTag','',@ischar);
addParameter(p,'dataS',[],@(x)validateattributes(x,{'struct'},{'nonempty'}));
addParameter(p,'dir_trainFile','',@(x)exist(x,'file'));
addParameter(p,'dir_dataFile','',@(x)exist(x,'file'));
addParameter(p,'gen_cond','',@(x)validateattributes(x,{'char'},{'nrows',1}));
addParameter(p,'gen_time','',@(x)any(validatestring(x,validGenTimes)));
addParameter(p,'progrMonitor',true,@islogical);
addParameter(p,'saveFile',true,@islogical);
addParameter(p,'trainS',[],@(x)validateattributes(x,{'struct'},{'nonempty'}));

% Parsing inputs.
parse(p,I);

% Assigning inputs to variables.
customTag = p.Results.customTag;
dataS = p.Results.dataS;
dirTrainFile = p.Results.dir_trainFile;
dirDataFile = p.Results.dir_dataFile;
genCond = p.Results.gen_cond;
genTimeStr = p.Results.gen_time;
progrMonitor = p.Results.progrMonitor;
saveFile = p.Results.saveFile;
trainS = p.Results.trainS;

% Custom checkings after parsing
requiredVars = {'dir_dataFile','gen_cond','gen_time'};
if any(ismember(requiredVars,p.UsingDefaults))
    error('mvpa:svmGeneralize:missingInput',...
        'All required parameters must be specified!');
end

[p,I] = deal([]); %#ok

%% Checking... 
% matlab version
if verLessThan('matlab', '8.3.0.532')
    error('mvpa:svmGeneralize:notSupportedMatlabVersion',...
        'MATLAB 2014a or newer is required! ');
end

% parallel pool 
currPool = gcp('nocreate');
if isempty(currPool)
    error('mvpa:svmGeneralize:missingParallelPool',...
        'Please open a parallel pool in order to run the function in parallel mode!');
end

% libsvm 
if isempty(regexp(path,'libsvm','once'))
    error('mvpa:svmGeneralize:missingToolbox','libsvm not found in path!');
end 

% parfor_progress
if progrMonitor
    if isempty(regexp(path,'parfor_progress','once'))
        error('mvpa:svmGeneralize:missingToolbox','parfor_progress not found in path!');
    end
end

%% Starting the timer
cStartFun = clock;
ds = datestr(now);
fprintf('\nsvmGeneralize%s%s\n%s\n',repmat(' ',1,67-length(ds)),ds,repmat('-',1,80));

%% Loading and checking data files...
if ~isempty(trainS)
    trDataset = trainS;
    trainS = []; %#ok
    dirTrainFile = trDataset.info.tr_filePath;
elseif ~isempty(dirTrainFile)
    trDataset = matfile(dirTrainFile);
else
    error('mvpa:svmGeneralize:missingInput',...
          'The training dataset is not specified!');
end

% Defining importand variables
trInfo = trDataset.info;
cvScheme = trInfo.cv_scheme;
scMethod = trInfo.sc_method;
scParam = trDataset.tr_scParam;
[~,name,ext] = fileparts(dirDataFile);
genDataFile = [name,ext];
if isfield(trInfo,'tr_data_file')
    trDataFile = trInfo.tr_data_file;
else
    % If the training data file is unknown we assume it is the same as the
    % generalization as this is the most likely case and also this is for
    % backward compatibility
    trDataFile = genDataFile;
end
trIsExample = trInfo.tr_isExample;
trMethod = trInfo.tr_method;
trTimePoints = trInfo.tr_timePoints;
trLabel = trInfo.tr_label;
trModels = trDataset.tr_models;
trGroupings = trDataset.tr_groupings;
% Number of training time points
nTrTimePoints = size(trTimePoints,2);
% Determining generalization label
genLabel = trLabel;
% Generalization time points
genTimePoints = trTimePoints;
% Cross-validation specific variables
if strcmp(cvScheme,'kf')
    nKFrep = trInfo.nKFrep;
    k = trInfo.k;
else
    nKFrep = [];
    k = [];
end
if regexp(trMethod,'-avg$','once')
    trAvgGrouping = trInfo.tr_avgGrouping;
    nRepsAvg = trInfo.tr_nRepsAvg;
else
    nRepsAvg = 1;
end

% Loading dataset to generalize to
if isempty(dataS)
    dataS = load(dirDataFile);
end
% dataM = matfile(dirDataFile);
dataMisc = dataS.misc;

[genPredlabels,genAccuracies,genGroupings] = deal(cell(nRepsAvg,1));

for iReps = 1:nRepsAvg
    % Selecting the appropriate set of examples and features for generalizing
    % genExamples_feat is a 3D matrix with 
    % size(1) = nFeatures
    % size(2) = nTimePoints
    % size(3) = nExamples
    % Feature scaling is also done here. 
    if ismember(trMethod,{'sample-wise-tf','sample-wise-tf-avg'})
        if regexp(trMethod,'-avg$','once')
            [genExamples_info,genExamples_feat,genIsExample] = ...
                mvpa.selectdata(dataS,genCond,trMethod,genLabel,genTimePoints,...
                                'tfFreq',trInfo.tf_freq,'tfType',trInfo.tf_type,...
                                'progrMonitor',progrMonitor,...
                                'trialGrouping',trAvgGrouping{iReps});
        else
            [genExamples_info,genExamples_feat,genIsExample] = ...
                mvpa.selectdata(dataS,genCond,trMethod,genLabel,genTimePoints,...
                                'tfFreq',trInfo.tf_freq,'tfType',trInfo.tf_type,...
                                'progrMonitor',progrMonitor);
        end
    else
        if regexp(trMethod,'-avg$','once')
            [genExamples_info,genExamples_feat,genIsExample] = ...
                mvpa.selectdata(dataS,genCond,trMethod,genLabel,genTimePoints,...
                                'progrMonitor',progrMonitor,...
                                'trialGrouping',trAvgGrouping{iReps});
        else
            [genExamples_info,genExamples_feat,genIsExample] = ...
                mvpa.selectdata(dataS,genCond,trMethod,genLabel,genTimePoints,...
                                'progrMonitor',progrMonitor);
        end
        
    end

    % Determining if the training condition and the generalization condition
    % are identical. 
    if ~strcmp(trDataFile,genDataFile)
        % If the training and generalization data files are different
        tr_genCondMatch = 'no';
    elseif all(~ismember(find(genIsExample),find(trIsExample)))
        % If the training and generalization example sets are disjoint
        tr_genCondMatch = 'no';
    elseif any(~ismember(find(trIsExample),find(genIsExample)))
        % If the training and generalization set is not disjoint but not all
        % training examples are present in the generalization example set. 
        tr_genCondMatch = 'partial';
        % error('mvpa:svmGeneralize:datasetMismatch',...
        % ['All training examples should be present among the generalization ',...
        % 'examples if they are not completely different!']);
    elseif all(ismember(find(genIsExample),find(trIsExample)))
        % If the generalization example set matches the training example set
        % completely
        tr_genCondMatch = 'full';
    else
        % If the generalization example set contains examples which are not
        % part of the training example set
        tr_genCondMatch = 'partial';
    end

    % Selecting generalization label
    if ~ismember(genLabel,genExamples_info.Properties.VariableNames)
        error('mvpa:svmGeneralize:invalidVariableValue',...
              'The specified generalization label is not part of the dataset!');
    else
        genExamples_lab = genExamples_info.(genLabel);
    end
    if any(isnan(genExamples_lab))
        error('mvpa:svmGeneralize:invalidVariableValue',...
              'Generalization labels can''t be NaNs! ');
    end
    sizeGenExamples_feat = size(genExamples_feat);

    % Determining the number of generalization time points. 
    if strcmp(genTimeStr,'tr')
        nGenTimePoints = 1;
    elseif strcmp(genTimeStr,'tr_x_tr')
        nGenTimePoints = nTrTimePoints;
    end
    
    trGroupingsActRep = trGroupings(:,:,iReps);
    
    % Grouping examples for generalization
    % For examples which are included in all folds (generalization to a
    % different condition) I use a label which is for sure not used as 
    % a label for in any other foldf (i.e. all fold labels are positive 
    % integers)
    commExmpLabel = -1;

    switch tr_genCondMatch
      case 'full'
        % The grouping is the same as for training
        genGroupingsActRep = trGroupingsActRep;
        
      case 'partial'
        % Finding training and generalization example IDs
        if regexp(trMethod,'-avg$','once')
            trIDs = unique(trAvgGrouping{iReps}(trIsExample & ~isnan(trAvgGrouping{iReps})));
            genIDs = unique(trAvgGrouping{iReps}(genIsExample & ~isnan(trAvgGrouping{iReps})));
        else
            trIDs = find(trIsExample);
            genIDs = find(genIsExample);
        end
        % Pre-allocating array for generalization groupings
        genGroupingsActRep = NaN(size(genExamples_info,1),size(trGroupingsActRep,2));
        genGroupingsActRep(ismember(genIDs,trIDs),:) = trGroupingsActRep(ismember(trIDs,genIDs),:);
        
        % Finding examples not present in the training example set
        if any(~ismember(genIDs,trIDs))
            % Assigning groups to the above examples according to the cv scheme
            if strcmp(cvScheme,'kf')
                misc.k = k;
                misc.nKFrep = nKFrep;
            end
            misc.nTimePoints = nTrTimePoints;
            tempGroupings = mvpa.assigngroupings(cvScheme,...
                                                 genExamples_info(~ismember(genIDs,trIDs),:),misc);
            genGroupingsActRep(~ismember(genIDs,trIDs),:) = tempGroupings;
        end
        
      case 'no'
        % % Assigning groups to the new examples according to the cv scheme
        % if strcmp(cvScheme,'kf')
        %     misc.k = k;
        %     misc.nKFrep = nKFrep;
        % end
        % misc.nTimePoints = nTrTimePoints;
        % genGroupingsActRep = mvpa.assigngroupings(cvScheme,genExamples_info,misc);
        genGroupingsActRep = ones(size(genExamples_info,1),1)*commExmpLabel;
    end

    % The maximum number of test examples in the cross-validation
    % scheme.
    nPredLabels = maxnival(genGroupingsActRep(:,1));
    
    % Selecting models for the actual repetition
    trModelsActRep = trModels(:,:,:,:,iReps);
    
    % Number of models per sample and in total. 
    nModelsPerTimePoint = size(trModelsActRep,3);
    
    scParamActRep = scParam(:,:,:,:,iReps);
    
    % Pre-allocating arrays for data collection. 
    genPredLabelsActRep = NaN(nPredLabels,nModelsPerTimePoint,nTrTimePoints,nGenTimePoints);
    genAccuraciesActRep = NaN(3,nModelsPerTimePoint,nTrTimePoints,nGenTimePoints);

    % Initializing progress monitor
    cStartGen = clock;
    fprintf('Generalizing...\n');
    if progrMonitor
        parfor_progress(nTrTimePoints*nGenTimePoints);
    end

    % Iterating through training samples
    parfor iTrTimePoint = 1:nTrTimePoints
        
        % Selecting the appropriate sample of models and scalin parameters. 
        actTrTimePointModels = trModelsActRep(:,:,:,iTrTimePoint);
        actTrTimePointScParams = scParamActRep(:,:,:,iTrTimePoint);
        
        % Pre allocating arrays for collecting data for a given sample. 
        actTrTimePointGenPredlab = NaN(nPredLabels,nModelsPerTimePoint,1,nGenTimePoints);
        actTrTimePointGenAcc = NaN(3,nModelsPerTimePoint,1,nGenTimePoints);
        
        for iGenTimePoint = 1:nGenTimePoints
            
            % Extracting the actual sample. This yields a matrix with
            % number of rows = number of examples,
            % number of columns = number of features.
            if strcmp(genTimeStr,'tr')
                actGenTimePointFeats = squeeze(genExamples_feat(:,iTrTimePoint,:))';
            else
                actGenTimePointFeats = squeeze(genExamples_feat(:,iGenTimePoint,:))';  %#ok
            end
            
            if any(size(actGenTimePointFeats) ~= [sizeGenExamples_feat(3),sizeGenExamples_feat(1)]) %#ok
                error('mvpa:svmGeneralize:datasetMismatch','Incorrect feature matrix sizes.');
            end
            
            % We do the generalization according to the cross-validation
            % scheme.
            
            % Preallocating arrays for test labels and features
            teLabs = NaN(nPredLabels,nModelsPerTimePoint);
            teFeats = NaN(nPredLabels,size(actGenTimePointFeats,2),nModelsPerTimePoint);
            
            if strcmp(cvScheme,'kf')
                for iKFrep = 1:nKFrep
                    actGrouping = genGroupingsActRep(:,iKFrep); %#ok
                    for iFold = 1:k
                        iModel = ((iKFrep-1)*k)+iFold;
                        actFoldExamples = ismember(actGrouping,[iFold,commExmpLabel]);
                        actNumExamples = sum(actFoldExamples);
                        teLabs(1:actNumExamples,iModel) = genExamples_lab(actFoldExamples); %#ok
                        teFeats(1:actNumExamples,:,iModel) = ...
                            mvpa.scalefeatures(actGenTimePointFeats(actFoldExamples,:),...
                                               scMethod,actTrTimePointScParams(:,:,iModel));
                    end
                end
            elseif strcmp(cvScheme,'loro')
                actGrouping = genGroupingsActRep;
                for iModel = 1:nModelsPerTimePoint
                    actModelExamples = ismember(actGrouping,[iModel,commExmpLabel]);
                    actNumExamples = sum(actModelExamples);
                    teLabs(1:actNumExamples,iModel) = genExamples_lab(actModelExamples);
                    teFeats(1:actNumExamples,:,iModel) = ...
                        mvpa.scalefeatures(actGenTimePointFeats(actModelExamples,:),...
                                           scMethod,actTrTimePointScParams(:,:,iModel));
                end
            end
            
            for iModel = 1:nModelsPerTimePoint
                
                actTeLabs = teLabs(:,iModel);
                actTeLabs = actTeLabs(~isnan(actTeLabs));
                actTeFeats = teFeats(:,:,iModel);
                actTeFeats = actTeFeats(~isnan(actTeLabs),:);
                
                [plabs,acc,~] = svmpredict(actTeLabs,actTeFeats,...
                                           mvpa.mat2mdl(actTrTimePointModels(:,:,iModel)));
                
                % Padding with NaNs if necesary
                if size(plabs,1) < nPredLabels
                    plabs = [plabs;NaN(nPredLabels-size(plabs,1),1)];
                end
                
                actTrTimePointGenPredlab(:,iModel,1,iGenTimePoint) = plabs;
                actTrTimePointGenAcc(:,iModel,1,iGenTimePoint) = acc;
                
            end
            
            % Advancing Progress monitor
            if progrMonitor
                parfor_progress;
            end
        end
        
        genPredLabelsActRep(:,:,iTrTimePoint,:) = actTrTimePointGenPredlab;
        genAccuraciesActRep(:,:,iTrTimePoint,:) = actTrTimePointGenAcc;
        
    end
    
    genPredlabels{iReps} = genPredLabelsActRep;
    genAccuracies{iReps} = genAccuraciesActRep;
    genGroupings{iReps} = genGroupingsActRep;
    
    % Finalizing progress monitor.
    if progrMonitor
        parfor_progress(0);
    end
    fprintf('Generalization elapsed time (days hours:minutes:seconds) %s \n\n',...
            datestr(etime(clock,cStartGen)/86400,'dd HH:MM:SS'));
end

% Saving the output into the struct. 
% Updating the info field
trInfo.gen_cond = genCond;
trInfo.gen_elapsedTime = etime(clock,cStartFun);
trInfo.gen_isExample = genIsExample;
trInfo.gen_label = genLabel;
trInfo.gen_time = genTimeStr;
[~,name,ext] = fileparts(dirDataFile);
trInfo.gen_data_file = [name,ext];
trInfo.gen_data_fs = dataMisc.fs;
% Saving computed data
genM = struct();
genM.gen_accuracies = mean(cat(5,genAccuracies{:}),5);
genM.gen_examples = genExamples_info;
genM.gen_groupings = cat(3,genGroupings{:});
if nRepsAvg > 1
    genM.gen_predlabels = catpad(5,genPredlabels{:});
else
    genM.gen_predlabels = genPredlabels{:};
end
genPredlabels = []; %#ok
genM.tr_accuracies = trDataset.tr_accuracies;
genM.tr_examples = trDataset.tr_examples;

% Closing the training file. 
trDataset = []; %#ok

% Ordering the sturcture fields. 
trInfo = orderfields(trInfo);
genM.info = trInfo; %#ok<*STRNU>

% Substring for generalization label
genLabelStr = genLabel;

% Generating output file name. 
if strcmp(genTimeStr,'tr')
    genTimeStr = 't';
elseif strcmp(genTimeStr,'tr_x_tr')
    genTimeStr = 'txt';
end

dateString = datestr(now,'yymmddHHMMSS');
if ~isempty(customTag), customTag = [customTag,'-']; end
fileNameAppend = sprintf('gen-%s%s_%s_%s_%s',customTag,genCond,genLabelStr,genTimeStr,...
    dateString);
extStart = strfind(dirTrainFile,'.mat');
filePath = [dirTrainFile(1:extStart-(length(dateString)+1)),...
    fileNameAppend,dirTrainFile(extStart:end)];
% Saving filename in the info
[~,fName,ext] = fileparts(filePath);
genM.info.gen_filename = [fName,ext];

% Generating the outputs. 
if saveFile
    save(filePath,'-struct','genM','-v7.3');
    varargout{1} = mvpares(filePath);
    varargout{2} = filePath;
else
    varargout{1} = mvpares(genM);
    varargout{2} = genM;
end


end