function prepdata(I)
% Prepares EEG data for mvpa. 
%
% INPUT: strucure I with the following fields   
%   dir_analysis:
%   dir_condDefFile:
%   dir_preproc:
%   
%   subID: 
%   tr_method: sample-wise, sample-wise-sm,...


%% Parsing input
p = inputParser;

requiredVars = {'dir_analysis','dir_preproc','subID','tr_method'};
validTrMethods = {'sample-wise','sample-wise-avg',...
    'sample-wise-bp','sample-wise-tf','sample-wise-tf-avg'};

addParameter(p,'dir_analysis','',@(x)exist(x,'dir'));
addParameter(p,'dir_preproc','',@(x)exist(x,'dir'));
addParameter(p,'subID','',@ischar);
addParameter(p,'tr_method','',@(x)any(validatestring(x,validTrMethods)));
addParameter(p,'euclNorm',true,@islogical);

parse(p,I);

if any(ismember(requiredVars,p.UsingDefaults))
    error('mvpa:prepdata:missingInput',...
        'All required parameters must be specified!');
end

dirAnalysis = p.Results.dir_analysis;
dirPreproc = p.Results.dir_preproc;
subID = p.Results.subID;
trMethod = p.Results.tr_method;
euclNorm = p.Results.euclNorm;

%% Checking setup
[~,setupID] = system('hostname');
setupID = regexp(setupID,'[a-zA-Z_0-9-]*','match');
setupID = setupID{1};

if ~isempty(regexp(setupID,'^bb.+','once'))
    isServer = true;
else
    isServer = false;
end

%%
% Load condition definition file
condDef = load(fullfile(get_path('project'),'experiment','trigger_LUT.mat'));
condDef = dataset2table(condDef.trigger_LUT);
if ~isa(condDef.blocktype,'categorical')
    condDef.blocktype = categorical(condDef.blocktype);
end

%% Extracting data from the available files
if strcmp(trMethod,'sample-wise')
    
    if euclNorm
        dataFname = [subID,'_','sw','_','data.mat'];
    else
        dataFname = [subID,'_','sw','_noEucl_','data.mat'];
    end
    eegFileIDstr = ['fteeg_MVPA_',subID,'_'];
    
    [feat,info,misc] = mvpa.extractdata(dirPreproc,eegFileIDstr,trMethod); %#ok
    
    % Multivariate noise normalization
    temp = permute(feat,[3,1,2]);
    aloc = unique(info.aloc);
    X = NaN(size(info,1),numel(aloc));
    for i = 1:numel(aloc)
        X(:,i) = ismember(info.aloc,aloc(i));
    end
    % X = [ones(size(info,1),1),info.aloc];
    sigma = mvpa.compCovEEG(temp,X);
    
    % Averaging the error covariance matrices across timepoints
    % sigma = mean(sigma,3);
    
    temp = mvpa.noiseNormalizeEEG(temp,sigma);
    feat = permute(temp,[2,3,1]);
    
    if euclNorm
        % Dividing the topographies by their Eucledian norm
        for i = 1:size(feat,2)
            for j = 1:size(feat,3)
                feat(:,i,j) = feat(:,i,j)./norm(feat(:,i,j));
            end
        end
    end
    
    % Smoothing the time courses of each trial
    % Window size of moving average, since the sampling rate is 1000 Hz,
    % sample numbers directly correspond to ms.
    % mavgwin = 20;
    % % Kernel for smoothing trials with the given time wintow
    % kern = ones(1,mavgwin)./mavgwin;
    
    % % if there is no parallel pool running, open one.
    % if isServer
    %     currPool = gcp('nocreate');
    %     if isempty(currPool)
    %         parpool('local',16);
    %     end
    % else
    %     currPool = gcp('nocreate');
    %     if isempty(currPool)
    %         parpool('local');
    %     end
    % end
    
    % for i = 1:size(feat,3)
        
    %     feattmp = feat(:,:,i);
    %     featsmtmp = feat(:,:,i);
        
    %     parfor j = 1:size(feattmp,1)
    %         featsmtmp(j,:) = conv(feattmp(j,:),kern,'same');
    %     end
        
    %     feat(:,:,i) = featsmtmp;
    % end
    
    % Saving the data
    save(fullfile(dirAnalysis,dataFname),'condDef','feat','info','misc','-v7.3');
    
elseif any(strcmp(trMethod,{'sample-wise-avg'}))
    % The averaging takes place within the training function, not
    % prior to training anymore, so the non averaged data are just
    % copied over
    if euclNorm
        sourceFname = [subID,'_','sw','_','data.mat'];
        destFname = [subID,'_','sw-avg','_','data.mat'];
    else
        sourceFname = [subID,'_','sw','_noEucl_','data.mat'];
        destFname = [subID,'_','sw-avg','_noEucl_','data.mat'];
    end
    
    copyfile(fullfile(dirPreproc,sourceFname),fullfile(dirAnalysis,destFname));
    
    % M = matfile(fullfile(dirPreproc,sourceFname));
    % info = M.info;
    % feat = M.feat;
    % misc = M.misc;
    
    % [info,infoAvg,featAvg,seed] = mvpa.averagetrials(info,feat,condDef,'rand',16);
    % % Preparing variables for saving
    % misc.info_source = info;
    % misc.seed = seed;
    % info = infoAvg;
    % feat = featAvg;
    % % Saving the data
    % save(fullfile(dirAnalysis,destFname),'condDef','feat','info','misc','-v7.3');
    
elseif strcmp(trMethod,'sample-wise-bp')
    
    destFname = [subID,'_','sw-bp','_','data.mat'];
    
    freqStr = {'d','t','a','b','gl','gh'};
    
    % Loading the first part of the dataset to determine the size. 
    eegFileIDstr = ['fteeg_MVPA-bp-',freqStr{1},'_'];
    [temp,~,~] = mvpa.extractdata(dirPreproc,eegFileIDstr,trMethod);
    tempSize = size(temp);
    temp = [];
    % Pre-allocating and saving the final array
    feat = NaN(tempSize(1),size(freqStr,2),tempSize(2),tempSize(3));
    save(fullfile(dirAnalysis,destFname),'feat','-v7.3');
    feat = [];
    % Opening the final array in writable mode
    M = matfile(fullfile(dirAnalysis,destFname),'Writable',true);
    % Pre-allocating array for misc information
    misc = cell(size(freqStr));
    for j = 1:numel(freqStr)
        eegFileIDstr = ['fteeg_MVPA-bp-',freqStr{j},'_'];
        [temp,info,misc{j}] = mvpa.extractdata(dirPreproc,eegFileIDstr,trMethod); %#ok
        % Saving the actual part of the array to disk
        M.feat(:,j,:,:) = reshape(temp,size(temp,1),1,size(temp,2),size(temp,3));
        temp = [];
    end
    M = [];
    % feat will be 4 dimensional:
    % 1. number of channels
    % 2. number of frequencies
    % 3. number of time samples
    % 4. number of examples
    temp = [misc{:}];
    misc = temp(1);
    misc.frequencies = {temp.frequencies};
    misc.freqStr = freqStr;
    
    % Appending the remaining variables
    save(fullfile(dirAnalysis,destFname),'condDef','info','misc','-append');
    
elseif strcmp(trMethod,'sample-wise-tf')
    
    destFname = [subID,'_','sw-tf','_','data.mat'];
    % Instantaneous power files
    eegFileIDstr = ['fteeg_MVPA_TF_',subID,'_'];
    [feat_p,info,misc] = mvpa.extractdata(dirPreproc,eegFileIDstr,trMethod); %#ok
    % Saving the first part of the file
    save(fullfile(dirAnalysis,destFname),'condDef','feat_p','info','misc','-v7.3');
    % Cleanin up memory from the big variable(s)
    feat_p = [];
    % Instantaneous phase files
    % eegFileIDstr = 'ptph_aeaMspmeeg_MVPA_TF';
    % [feat_ph,~,~] = mvpa.extractdata(dirPreproc,eegFileIDstr,trMethod); %#ok
    % save(fullfile(dirAnalysis,destFname),'feat_ph','-append');
    
    % feat_p and feat_ph will be 4 dimensional: 
    % 1. number of channels
    % 2. number of frequencies
    % 3. number of time samples
    % 4. number of examples
elseif strcmp(trMethod,'sample-wise-tf-avg')
    
    sourceFname = [subID,'_','sw-tf','_','data.mat'];
    destFname = [subID,'_','sw-tf-avg','_','data.mat'];

    M = matfile(fullfile(dirPreproc,sourceFname));
    info = M.info;
    feat_p = M.feat_p;
    misc = M.misc;
    
    [info,infoAvg,featAvg,seed] = mvpa.averagetrials(info,feat_p,condDef,'rand');
    % Preparing variables for saving
    misc.info_source = info;
    misc.seed = seed;
    info = infoAvg;
    feat_p = featAvg;
    % Saving the data
    save(fullfile(dirAnalysis,destFname),'condDef','feat_p','info','misc','-v7.3');
        
end

end


function info = assignnewblocks(info,condDef)

blockTypes = unique(info.blocktype);
blocks = unique(info(:,{'session','run','block','blocktype'}),'rows');
blocks.blockID = (1:size(blocks))';

actConds = condDef.condition(ismember(condDef.blocktype,blockTypes));
nConds = size(actConds,1);
info.blockIDnew = NaN(size(info,1),1);
for j = 1:nConds
    isExample = info.condition == actConds(j) & ~info.catch_trial;
    if all(~isExample)
        error('mvpa:prepdata:missingExample',...
                'No example was found for condition %d, ',...
                 actConds(j));
    end
    actBlocks = blocks(blocks.blocktype == condDef.blocktype(...
        condDef.condition == actConds(j)),:);
    nBlocks = size(actBlocks,1);
    blockIDnew = mod(randperm(sum(isExample))',nBlocks);
    blockIDnew(blockIDnew == 0) = nBlocks;
    for i = 1:size(blockIDnew), blockIDnew(i) = actBlocks.blockID(blockIDnew(i)); end
    info.blockIDnew(isExample) = blockIDnew;
end

end
