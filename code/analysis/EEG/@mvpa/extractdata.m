function [feat,info,misc] = extractdata(dataDir,idStr,tr_method)
% Extracts the necessary data for MVPA for a given subject. 
%
% Input: 
%   dataDir: the full path of the data folder to work on. 
%   idStr: The identification string which determines which eeg files are
%       to be used. This is the first part of the name of the file ahead of
%       the first '_' character. 
%   tr_method: training method
% Output: 
%   info: table containing all information necessary from the examples
%       for MVPA. The size is the number of examples x 5 variables
%       Variables (in order): 
%           session
%           condition
%           task
%           locA
%           locV
%   feat: L x M x N 3D or L x M x N x O 4D array containing the features, 
%       if 3D: 
%           L is the number of channels(features)
%           M is the number of time samples
%           N is the number of examples 
%       if 4D: 
%           L is the number of channels(features)
%           M is the number of frequencies
%           N is the number of time samples
%           O is the number of examples 
%
% Details: 
%   Trials marked as 'bad' during the preprocessing are ignored. 
%   size(data.info,1) == N should be true! 
%
%% Parsing input
p = inputParser;

validTrMethods = {'sample-wise','sample-wise-bp','sample-wise-tf'};

addRequired(p,'dataDir',@ischar);
addRequired(p,'idStr',@ischar);
addRequired(p,'tr_method',@(x)any(validatestring(x,validTrMethods)));

parse(p,dataDir,idStr,tr_method);

dataDir = p.Results.dataDir;
idStr = p.Results.idStr;
tr_method = p.Results.tr_method;

%%
saveDF = cd(dataDir);

% list matching files
fileNameList = cellstr(ls([idStr '*.mat']));
if isempty(fileNameList)
    error('mvpa:extractdata:invalidInput',...
        'No eeg files were found for feature extraction! ');
end

filePathList = fullfile(repmat({pwd},size(fileNameList)),fileNameList);

cd(saveDF);

[infos,infosAll,feats] = deal(cell(1,size(fileNameList,1)));

for i = 1:size(filePathList,1)
    
    % Load fieldtrip data file
    ftData = load(filePathList{i});
    fieldNames = fieldnames(ftData);
    if numel(fieldNames) == 1
        ftData = ftData.(fieldNames{1});
    else
        error('mvpa:extractdata:invalidInput',...
            'Data files must a single fieldtrip data structure');
    end
    
    % get trial info
    if ~isfield(ftData,'trialinfo')
        error('mvpa:extractdata:missingField',...
            'FieldTrip data must contain the trialinfo field');
    end
        
    tempInfo = struct2table(cell2mat(ftData.trialinfo));
    tempInfo.blocktype = categorical(tempInfo.blocktype);
    tempInfo.adaptdir = categorical(tempInfo.adaptdir);
    tempInfo.hand = categorical(tempInfo.hand);
    % logical vector for good trials
    isGoodTrial = true(size(tempInfo,1),1);
    isGoodTrial(tempInfo.badtrials ~= 0) = false;
    
    
    % getting the features
    
    if strcmp(ft_datatype(ftData),'raw')
        tempFeats = cat(3,ftData.trial{:});
        tempFeats = tempFeats(:,:,isGoodTrial);
    else
        [lia,locb] = ismember({'chan','freq','time','rpt'}, ...
                              regexp(ftData.dimord,'_','split'));
        if all(lia)
            tempFeats = ftData.powspctrm;
            tempFeats = permute(tempFeats,locb);
            tempFeats = tempFeats(:,:,:,isGoodTrial);
        else
            error('mvpa:extractdata:missingDimension',...
                  'FieldTrip data does not contain a required dimension');
        end
    end
    
    infosAll{i} = tempInfo;
    tempInfo.badtrials = [];
    infos{i} = tempInfo(isGoodTrial,:);
    feats{i} = tempFeats;
    
    % Free up memory
    tempInfo = [];
    tempFeats = [];
end

if isfield(ftData,'fsample')
    misc.fs = ftData.fsample;
elseif isfield(ftData,'time')
    if strcmp(ft_datatype(ftData),'raw')
        misc.fs = round(1/(ftData.time{1,1}(2)-ftData.time{1,1}(1)));
    else
        misc.fs = round(1/(ftData.time(2)-ftData.time(1)));
    end
else
    error('mvpa:extractdata:missingField',...
        'FieldTrip data must contain the time field');
end
if strcmp(ft_datatype(ftData),'raw')
    misc.timeOnset = ftData.time{1,1}(1);
else
    misc.timeOnset = ftData.time(1);
end
misc.featureLabels = ftData.label;
misc.infoAll = cat(1,infosAll{:});

if strcmp(tr_method,'sample-wise-tf')
    if isfield(ftData,'freq')
        misc.frequencies = ftData.freq;
    else
        error('mvpa:extractdata:missingField',...
              ['FieldTrip data must contain the cfg.freq ' ...
               'field.']);
    end
elseif strcmp(tr_method,'sample-wise-bp')
    cfg = ftData.cfg;
    while 1
        if isfield(cfg,'bpfreq')
            misc.frequencies = cfg.bpfreq;
            break;
        elseif isfield(cfg,'previous')
            cfg = cfg.previous;
        else
            error('mvpa:extractdata:missingField',...
                'FieldTrip data must contain the cfg.bpfreq field.');
        end
    end
end
featSize = size(feats{1,1});

info = cat(1,infos{:});
[infos,infosAll,ftData] = deal([]);

if size(featSize,2) > 3
    feat = cat(4,feats{:});
    feats = [];
    if size(info,1) ~= size(feat,4)
        error('Inconsistent number of examples! ');
    end
else
    feat = cat(3,feats{:});
    feats = [];
    if size(info,1) ~= size(feat,3)
        error('Inconsistent number of examples! ');
    end
end

end
