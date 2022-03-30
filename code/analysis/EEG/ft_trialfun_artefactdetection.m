function [trl,event] = ft_trialfun_artefactdetection(cfg)

% Finding the actual trials in the data based on the specified event values
hdr = ft_read_header(cfg.headerfile);
event = ft_read_event(cfg.headerfile);
evValues = {event.value}';
evTypes = {event.type}';
evSamples = [event.sample]';

isReqType = strcmp(evTypes,cfg.trialdef.eventtype);
evValues = evValues(isReqType);
evSamples = evSamples(isReqType);

% trigger values
trStim = cfg.trialdef.trigdef.trig_eeg(ismember(cfg.trialdef.trigdef.type,'stim'));
trBlockStart = cfg.trialdef.trigdef.trig_eeg(ismember(cfg.trialdef.trigdef.type,'blockstart'));

% Correcting for trigger delay with respect to visual onset
temp = evSamples(ismember(evValues,[trStim',trBlockStart']));
temp = temp+(cfg.trialdef.trig_visonset_corr*hdr.Fs);
evSamples(ismember(evValues,[trStim',trBlockStart'])) = temp;

% all event values must be strings
nonStrIdx = find(~cellfun(@isstr,evValues));
evValues(nonStrIdx) = repmat({''},numel(nonStrIdx),1);

% Checking the number of sessions
nBlocks = sum(ismember(evValues,trBlockStart));
if nBlocks == 0
    error('No session start triggers were found!');
elseif nBlocks ~= cfg.trialdef.fileSpec.nBlocksInFile
    error('Number of sessions does not match the expected!');
end
blockStartEvSamples = evSamples(ismember(evValues,trBlockStart));
stimOnsetEvSamples = evSamples(ismember(evValues,trStim));

% Removing blocks which are to be excluded
blockToExclude = cfg.trialdef.fileSpec.blocksToExclude;
if ~isnan(blockToExclude)
    for i = 1:numel(blockToExclude)
        % Find and remove trials which belong to the block(s) to be
        % removed. 
        if blockToExclude(i) ~= nBlocks
            startExcluded = blockStartEvSamples(blockToExclude(i));
            endExcluded = blockStartEvSamples(blockToExclude(i)+1);
            stimOnsetEvSamples(stimOnsetEvSamples >= startExcluded & stimOnsetEvSamples < endExcluded) = [];
        else
            startExcluded = blockStartEvSamples(blockToExclude(i));
            stimOnsetEvSamples(stimOnsetEvSamples >= startExcluded) = [];
        end
    end
    % Removing blocks to be excluded and updating number of blocks
    blockStartEvSamples(blockToExclude) = [];
    nBlocks = numel(blockStartEvSamples);
end

se = NaN(nBlocks,2);

% Finding session start and end samples. 
if nBlocks == 1
    
    se(1) = blockStartEvSamples;
    se(2) = stimOnsetEvSamples(end)+round(cfg.trialdef.poststim*hdr.Fs);
    
else
    
    se(:,1) = blockStartEvSamples;
    
    for i = 1:nBlocks
        
        % Find the last stimulus onset in the given session
        if i == nBlocks
            lastStimInBlockSample = stimOnsetEvSamples(end);
        else
            lastStimInBlockSample  = max(stimOnsetEvSamples(stimOnsetEvSamples < se(i+1,1)));
        end
        se(i,2) = lastStimInBlockSample+round(cfg.trialdef.poststim*hdr.Fs);
        
    end

end

% Generating fake trials of length ~cfg.trialdef.trllength just for the 
% efficient artefact reviewing.
trlLengthSampl = cfg.trialdef.faketrllength*hdr.Fs;

% Number of artificial trials for each session. 
nTrialArtfPerBlock = ceil((se(:,2)-se(:,1))/trlLengthSampl);

nTrialArtf = sum(nTrialArtfPerBlock);

trl = zeros(nTrialArtf,3);

iTrial = 1;
for iBlock = 1:nBlocks
    
    startIdx = se(iBlock,1);
    
    for iTrialArtfPerBlock = 1:nTrialArtfPerBlock(iBlock)
        
        % Marking the beginning and ending samples of trials for artefact
        % reviewing. The offset remains always zero.
        trl(iTrial,1) = startIdx;
        trl(iTrial,2) = startIdx+(trlLengthSampl-1);
        startIdx = startIdx+trlLengthSampl;
        iTrial = iTrial+1;
        
    end
end

end