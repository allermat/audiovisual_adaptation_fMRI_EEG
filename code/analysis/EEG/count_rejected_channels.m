% This script reads the pre-processed EEG data and extracts which channels
% were rejected and interpolated
clearvars;
subIDlist = {'149','336','340','345','346'};
badChannels = cell(size(subIDlist));
for iSub = 1:numel(subIDlist)
    subID = subIDlist{iSub};
    sourceDir = fullfile(get_path('project'),'data',subID,'EEG','preproc_data','COMM');
    listing = dir(fullfile(sourceDir,'fteeg_*.mat'));
    ftData = cellfun(@load,fullfile(sourceDir,{listing.name}'));
    % Some FT files have an elec field, some not, so can't use regular
    % concatenation
    temp = catpadstruct(ftData.ftDataReref);
    temp = [temp.cfg]';
    % Whenever there was no bad channel, the channel interpolation was not
    % done, so the cfg.previous field points to the filtering step. The
    % padded concatenation works well here as the badchannel field is only
    % populated whenever bad channels were rejected. 
    cfg = catpadstruct(temp.previous);
    badChannels{iSub} = {cfg.badchannel};
end

% Quantification per file
for i = numel(badChannels):-1:1
    badChanStats.max(i) = max(cellfun(@numel,badChannels{i}));
    badChanStats.min(i) = min(cellfun(@numel,badChannels{i}));
    badChanStats.mean(i) = mean(cellfun(@numel,badChannels{i}));
end

% Quantification per subject, pooled over files
cellfun(@(c) unique(cat(1,c{:})),badChannels,'UniformOutput',false)