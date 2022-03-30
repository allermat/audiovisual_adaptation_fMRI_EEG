function gr_analysis_ERP()
% Group level ERP analysis

saveDf = cd(fullfile(get_path('project'),'data'));
subjList = cellstr(ls);
subjList = regexp(subjList,'[0-9]{3}','match','once');
subjList = subjList(~ismember(subjList,{''}));
cd(saveDf);

ftFilesToAvg = cell(size(subjList));
matchStrTokens = 'fteeg_ERP_[0-9]{3}(.*).mat';

% Finding match strings for conditions
saveDf = cd(fullfile(get_path('project'),'data',subjList{1},'EEG','preproc_data','ERP'));
fileList = cellstr(ls);
cd(saveDf);
temp = regexp(fileList,matchStrTokens,'tokens');
temp = temp(~cellfun(@isempty,temp));
temp = [temp{:}]';
matchStrConds = [temp{:}]';

% Grand averages separately for aloc, pretest, posttest L, posttest R
for i = 1:size(matchStrConds,1)
    
    for j = 1:size(subjList,1)
        
        saveDf = cd(fullfile(get_path('project'),'data', ...
            subjList{j},'EEG','preproc_data','ERP'));
        fileList = cellstr(ls);
        matchID = ~cellfun(@isempty,regexp(fileList,[matchStrConds{i},'.mat']));
        
        if sum(matchID) == 0
            warning('No file, skipping this subject! ');
            cd(saveDf);
            continue;
        elseif sum(matchID) > 1
            warning('Multiple files, skipping this subject! ');
            cd(saveDf);
            continue;
        else
            ftFilesToAvg{j} = load(fileList{matchID});
            ftFilesToAvg{j} = ftFilesToAvg{j}.ftDataAvg;
        end
        
        cd(saveDf);
    end
    
    ftFilesToAvg = ftFilesToAvg(~cellfun(@isempty,ftFilesToAvg));
    
    % Taking grand average across subjects
    cfg = struct();
    ftDataGrAvg = ft_timelockgrandaverage(cfg,ftFilesToAvg{:});
    % Getting rid of the unnecessary previous nested configs
    ftDataGrAvg.cfg.previous = [];
    % Saving data
    fprintf('\n\nSaving data...\n\n');
    fileName = ['fteeg_gr_avg_ERP',matchStrConds{i},'.mat'];
    savePath = fullfile(get_path('project'),'data', ...
        'group','EEG','ERP',fileName);
    save(savePath,'ftDataGrAvg','-v7.3');
    
end

end

