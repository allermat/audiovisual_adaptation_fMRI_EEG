function [feat,info] = pcm_getData_EEG(timeWins,output,varargin) 

p = inputParser;

addRequired(p,'timeWins');
addRequired(p,'output');
addParameter(p,'blocktype','pre-postL-postR');
addParameter(p,'iseucnorm',false,@(x) validateattributes(x,{'logical'},{'scalar'}));

parse(p,timeWins,output,varargin{:});

timeWins = p.Results.timeWins;
output = p.Results.output;
blocktype = p.Results.blocktype;
iseucnorm = p.Results.iseucnorm;

subID = {'149','336','340','345','346'};
[temp,info] = deal(cell(numel(subID),1));
feat = cell(numel(subID),numel(timeWins));
trainMethod = 'sample-wise-avg';
timePoints = -100:5:500;
if strcmp(blocktype,'pre')
    condSel = blocktype;
    rdmLayout = blocktype;
elseif strcmp(blocktype,'pre-postL-postR')
    condSel = 'all'; 
    rdmLayout = blocktype;
end
for iSubj = 1:numel(subID)
    if iseucnorm
        dataFileName = [subID{iSubj},'_sw-avg_data.mat'];
    else
        dataFileName = [subID{iSubj},'_sw-avg_noEucl_data.mat'];
    end
    
    I = struct();
    I.dir_dataFile = fullfile(get_path('project'),'data',subID{iSubj},...
                              'EEG','PCM',trainMethod,dataFileName);
    I.condSel = condSel;
    I.rdmLayout = rdmLayout;
    I.tr_method = trainMethod;
    I.timePoints = num2cell(timePoints);
    I.output = output;
    [temp{iSubj},info{iSubj}] = pcm_extractdata(I);
    
    for iTimeWin = 1:numel(timeWins)
        isInTimeWin = ismember(timePoints,timeWins{iTimeWin});
        feat{iSubj,iTimeWin} = squeeze(mean(temp{iSubj}(:,isInTimeWin,:),2))';
    end
end

info = repmat(info,1,size(feat,2));

end