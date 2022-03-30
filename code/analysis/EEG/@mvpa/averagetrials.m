function [info,infoAvg,featAvg,seed] = averagetrials(info,feat,condDef,mode,varargin)
% Method for averaging trials of the same condition within each block (period)

% Parsing input
p = inputParser;

validModes = {'rand','run','preset'};

addRequired(p,'info',@(x) validateattributes(x,{'table'},{'nonempty'}));
addRequired(p,'feat',@(x) validateattributes(x,{'numeric'},{'nonempty'}));
addRequired(p,'condDef',@(x) validateattributes(x,{'table'},{'nonempty'}));
addRequired(p,'mode',@(x)any(validatestring(x,validModes)));
addOptional(p,'nTrialsToAvg',8,...
            @(x) validateattributes(x,{'numeric'},...
                                    {'scalar','integer','>',1}));
addOptional(p,'trialGrouping',[],...
            @(x) validateattributes(x,{'numeric'},...
                                    {'column','nrows',size(info,1)}));

parse(p,info,feat,condDef,mode,varargin{:});

info = p.Results.info;
feat = p.Results.feat;
condDef = p.Results.condDef;
mode = p.Results.mode;
nTrialsToAvg = p.Results.nTrialsToAvg;
trialGrouping = p.Results.trialGrouping;

if strcmp(mode,'preset') && isempty(trialGrouping)
    error('mvpa:averagetrials:missingInput',...
          'trialGrouping must be specified if the mode is ''preset''');
end

featSize = num2cell(size(feat));
dimToAvg = numel(featSize);

if ismember(mode,{'rand','preset'})
    if strcmp(mode,'rand')
        [info,seed] = grouptrials(info,condDef,nTrialsToAvg);
    else
        info.avgExampleID = trialGrouping;
        seed = [];
    end
    [avgExampleIDs,id] = unique(info.avgExampleID(~isnan(info.avgExampleID)));
    featAvg = NaN(featSize{1:dimToAvg-1},numel(avgExampleIDs));
    nTrialsAvgd = NaN(numel(avgExampleIDs),1);
    for j = 1:numel(avgExampleIDs)
        isExample = info.avgExampleID == avgExampleIDs(j);
        if dimToAvg == 3
            featAvg(:,:,j) = mean(feat(:,:,isExample),dimToAvg);
        elseif dimToAvg == 4
            featAvg(:,:,:,j) = mean(feat(:,:,:,isExample),dimToAvg);
        end
        nTrialsAvgd(j) = sum(isExample);
    end
    
    % Constructing info
    tempInfo = info(~isnan(info.avgExampleID),:);
    tempInfo = tempInfo(id,:);
    tempInfo.iExample = tempInfo.avgExampleID;
    tempInfo = tempInfo(:,[end,1:end-1]);
    tempInfo(:,{'session','run','block','avgExampleID','resp','RT','hand'}) = [];
    tempInfo.nTrialsAvgd = nTrialsAvgd;
    infoAvg = tempInfo;

else
    featAvg = []; %NaN(featSize{1:dimToAvg-1},nBlocks*nConds);
    infoAvg = table();
    blocks = unique(info(:,{'session','run','block'}),'rows');
    nBlocks = size(blocks,1);
    seed = [];
    for i = 1:nBlocks
        actBlockInfo = info(ismember(info(:,{'session','run','block'}),blocks(i,:)),:);
        actConds = condDef.condition(condDef.blocktype == actBlockInfo.blocktype(1));
        nConds = size(actConds,1);
        temp = NaN(featSize{1:dimToAvg-1},nConds);
        nTrialsAvgd = NaN(nConds,1);
        missingCond = false;
        for j = 1:nConds
            isExample = ismember(info(:,{'session','run','block'}),blocks(i,:)) & ...
                info.condition == actConds(j) & ~info.catch_trial;
            if all(~isExample)
                warning('mvpa:averagetrials:missingExample',...
                        ['No example was found for block %d condition %d, ' ...
                         'skipping this block!'],i,actConds(j));
                missingCond = true;
                break;
            end
            if dimToAvg == 3
                temp(:,:,j) = mean(feat(:,:,isExample),dimToAvg);
            elseif dimToAvg == 4
                temp(:,:,:,j) = mean(feat(:,:,:,isExample),dimToAvg);
            end
            nTrialsAvgd(j) = sum(isExample);
        end

        if ~missingCond
            % Getting rid of the catch trials
            actBlockInfo = actBlockInfo(~actBlockInfo.catch_trial,:);
            [~,ia] = unique(actBlockInfo(:,{'session','run','block','condition'}),'rows');
            tempInfo = actBlockInfo(ia,:);
            if any(tempInfo.condition ~= actConds)
                error('mvpa:averagetrials:conditionMismatch',...
                      'Inconsistent order of conditions! ');
            end
            tempInfo.nTrialsAvgd = nTrialsAvgd;
            infoAvg = cat(1,infoAvg,tempInfo);
            featAvg = cat(dimToAvg,featAvg,temp);
        end
    end
end
if ismember('blockIDnew',infoAvg.Properties.VariableNames)
    infoAvg.blockIDnew = [];
end
% Checking if the respective sizes match
if size(infoAvg,1) ~= size(featAvg,dimToAvg)
    error('mvpa:averagetrials:exampleNumberMismatch',...
        'Inconsistent number of examples in info and feat! ');
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
factors = fullfact([size(actConds,1),numel(adaptDir), ...
                    numel(catchLevels)]);
% factors = fullfact([size(actConds,1),numel(sessions),numel(catchLevels)]);

rng('shuffle');
seed = rng;

for iFact = 1:size(factors,1)
    % Making sure only trials with same response hand are
    % averaged 
    isExample = info.condition == actConds(factors(iFact,1)) & ...
        info.adaptdir == adaptDir(factors(iFact,2)) & ...
        info.catch_trial == catchLevels(factors(iFact,3));
    if all(~isExample)
        warning('mvpa:averagetirals:grouptrials:missingExample',...
              'No example was found for condition %d, ', ...
              actConds(factors(iFact,1)));
    end

    % Number of examples after averaging in this condition
    nExamplesOut = floor(sum(isExample)/nTrialsToAvg);
    if nExamplesOut < 1
        warning('mvpa:averagetrials:grouptrials:missingExample',...
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