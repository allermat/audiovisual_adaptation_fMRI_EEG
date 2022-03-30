function [trl,event] = ft_trialfun_eventlocked(cfg)

hdr = ft_read_header(cfg.headerfile);
event = ft_read_event(cfg.headerfile);
evValues = {event.value}';
evTypes = {event.type}';
evSamples = [event.sample]';

trigDef = cfg.trialdef.trigdef;

isReqType = strcmp(evTypes,cfg.trialdef.eventtype);
evValues = evValues(isReqType);
evSamples = evSamples(isReqType);

% Blocktypes
allBlockTypes = unique(trigDef.blocktype(ismember(trigDef.type,'stim')));

if ismember(cfg.trialdef.blocktype,'all')
    blocksToUse = allBlockTypes;
else
    blocksToUse = cfg.trialdef.blocktype;
end
if any(~ismember(blocksToUse,allBlockTypes))
    error('ft_trialfun_eventlocked:invalidBlockType',...
          'The specified blocktype is invalid');
end

% trigger values
trStim = trigDef.trig_eeg(...
    ismember(trigDef.type,'stim') & ...
    ismember(trigDef.blocktype,blocksToUse));
trBlockStart = ...
    trigDef.trig_eeg(ismember(trigDef.type,'blockstart'));
trRunStart = ...
    trigDef.trig_eeg(ismember(trigDef.type,'expstart'));
trRunStop = ...
    trigDef.trig_eeg(ismember(trigDef.type,'expend'));

% Correcting for trigger delay with respect to visual onset
temp = evSamples(ismember(evValues,[trStim',trBlockStart',trRunStart',trRunStop']));
temp = temp+(cfg.trialdef.trig_visonset_corr*hdr.Fs);
evSamples(ismember(evValues,[trStim',trBlockStart',trRunStart',trRunStop'])) = temp;

% all event values must be strings
nonStrIdx = find(~cellfun(@isstr,evValues));
evValues(nonStrIdx) = repmat({''},numel(nonStrIdx),1);

% Checking the number of blocks
nRuns = sum(ismember(evValues,trRunStart));
if nRuns == 0
    error('No run start triggers were found!');
elseif nRuns ~= cfg.trialdef.fileSpec.nRunsInFile
    error('Number of runs does not match the expected!');
end
% Finding various event samples
targEvSamples = evSamples(ismember(evValues,trStim));
targEvValues = evValues(ismember(evValues,trStim));
runStartEvSamples = evSamples(ismember(evValues,trRunStart));
% I intentionally don't use the runStopEvSamples because if the run
% was aborted (externally, not by quitting) then there is just a
% run start trigger, but no run stop trigger. 
% runStopEvSamples = evSamples(ismember(evValues,trRunStop));
blockStartEvSamples = evSamples(ismember(evValues,trBlockStart));

% Removing runs which are to be excluded
runToExclude = cfg.trialdef.fileSpec.runsToExclude;
if ~isnan(runToExclude)
    for i = 1:numel(runToExclude)
        % Find and remove trials which belong to the run(s) to be
        % removed. 
        if runToExclude(i) ~= nRuns
            startExcluded = runStartEvSamples(runToExclude(i));
            endExcluded = runStartEvSamples(runToExclude(i)+1);
            % Clearing stimulus and other event samples and values
            % within the run(s) to be excluded
            targEvValues(targEvSamples >= startExcluded & ...
                         targEvSamples < endExcluded) = [];
            targEvSamples(targEvSamples >= startExcluded & ...
                          targEvSamples < endExcluded) = [];
            blockStartEvSamples(blockStartEvSamples >= startExcluded & ...
                                blockStartEvSamples < endExcluded) = [];
        else
            startExcluded = runStartEvSamples(runToExclude(i));
            % Clearing stimulus and other event samples and values
            % within the run(s) to be excluded
            targEvValues(targEvSamples >= startExcluded) = [];
            targEvSamples(targEvSamples >= startExcluded) = [];
            blockStartEvSamples(blockStartEvSamples >= startExcluded) = [];
        end
    end
    % Removing runs to be excluded and updating number of runs
    runStartEvSamples(runToExclude) = [];
    % runStopEvSamples(runToExclude) = [];
end
% Grouping block start event samples by runs
blockStartEvSamplesInRuns = cell(size(runStartEvSamples));
for i = 1:size(runStartEvSamples)
    if i < size(runStartEvSamples)
        blockStartEvSamplesInRuns{i} = blockStartEvSamples(...
            blockStartEvSamples >= runStartEvSamples(i) & ...
            blockStartEvSamples < runStartEvSamples(i+1));
    else
        blockStartEvSamplesInRuns{i} = blockStartEvSamples(...
            blockStartEvSamples >= runStartEvSamples(i));
    end
end
nTrials = numel(targEvSamples);

[begSamples,endSamples] = deal(NaN(nTrials,1));
[run,block,iTrialInBlock,cond] = deal(NaN(nTrials,1));

% Initializing counters
actRunInFile = 0;
actBlockInRun = 0;
actTrialInBlock = 0;

for i = 1:nTrials
    
    % Find actual target event sample
    actTargEvSampl = targEvSamples(i);
    % Updating actual run
    if actRunInFile ~= find(actTargEvSampl > runStartEvSamples,1,'last')
        actRunInFile = find(actTargEvSampl > runStartEvSamples,1,'last');
        % Resetting block counter
        actBlockInRun = 0;
    end
    % Updating actual block and actual trial in block
    if actBlockInRun ~= find(actTargEvSampl > blockStartEvSamplesInRuns{actRunInFile},1,'last')
        actBlockInRun = find(actTargEvSampl > blockStartEvSamplesInRuns{actRunInFile},1,'last');
        % Resetting trial counter
        actTrialInBlock = 0;
    end
    % Updating trial counter
    actTrialInBlock = actTrialInBlock + 1;
    % Beginning and ending samples of trials
    begSamples(i) = actTargEvSampl-round(cfg.trialdef.prestim*hdr.Fs);
    endSamples(i) = actTargEvSampl+round(cfg.trialdef.poststim*hdr.Fs);
    % Run number
    run(i) = cfg.trialdef.fileSpec.runid(actRunInFile);
    % Block number
    block(i) = actBlockInRun;
    % Serial number of trial in the block
    iTrialInBlock(i) = actTrialInBlock;
    % Stimulus condition
    cond(i) = trigDef.condition(ismember(trigDef.trig_eeg,targEvValues{i}));
    
end

% Creating the offset array
offset = -round(cfg.trialdef.prestim*hdr.Fs)*ones(size(begSamples));
% Creating the session array
session = cfg.trialdef.fileSpec.session*ones(size(begSamples));
trl = [begSamples,endSamples,offset,session,run,block,iTrialInBlock,cond];

end
