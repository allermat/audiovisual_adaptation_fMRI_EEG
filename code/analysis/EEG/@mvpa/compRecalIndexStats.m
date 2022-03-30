function stats = compRecalIndexStats(indivData,time,varargin)
% Method to compute statistics for recalibration index.

% Parsing input
p = inputParser;
addRequired(p,'indivData',@(x) validateattributes(x,{'struct'},{'vector',...
    'nonempty'}));
addRequired(p,'time',@(x) validateattributes(x,{'numeric'},{'vector',...
    'increasing'}));
addOptional(p,'timeWin',{},@(x) validateattributes(x,{'cell'},{'nonempty'}));
parse(p,indivData,time,varargin{:});
indivData = p.Results.indivData;
time = p.Results.time;
timeWin = p.Results.timeWin;

if isrow(time), time = time'; end

if isempty(timeWin)
    isInTimeWin = {true(size(time))};
else
    fun = @(x) time >= x(1) & time <= x(2);
    isInTimeWin = cellfun(fun,timeWin,'UniformOutput',false);
    % Checking if the timewindows overlap
    temp = cat(2,isInTimeWin{:});
    if any(sum(temp,2) > 2)
        error('mvpa:compRecalIndexStats:invalidInput',...
              'The time windows of interest must not overlap!');
    end
end

% Number of subjects
nSubj = length(indivData);
% Analysis names.
anNames = {'RI','RIbin'};
anVarNames = {'RI','RIbin'};
tempStats = cell(size(anNames));
% Full factorial expanson of the factor levels always follows this patern 
ff = [1,2];
nCond = size(ff,2);
dataFieldNames = fieldnames(indivData);
indivData = squeeze(struct2cell(indivData));
for i = 1:numel(anNames)
    % Choosing variables for the actual analysis
    [~,actVarIdx] = ismember(anVarNames{i},dataFieldNames);
    temp = indivData(actVarIdx,:);
    data = cat(2,temp{:});
    [tempStats{i}.prob,tempStats{i}.stat,tempStats{i}.h,...
    tempStats{i}.prob_uc,tempStats{i}.h_uc] = deal(NaN(size(time)));
    for iTimeWin = 1:numel(isInTimeWin)
        % Computing statistics
        [pu,hu,pc,hc,st] = mvpa.bootstrpOneSampleTtest(data(isInTimeWin{iTimeWin},:),...
                                                     10000,'MCPsol','cluster',...
                                                     'clusterStat','maxsize');
        tempStats{i}.prob(isInTimeWin{iTimeWin}) = pc;
        tempStats{i}.h(isInTimeWin{iTimeWin}) = hc;
        tempStats{i}.stat(isInTimeWin{iTimeWin}) = st;
        tempStats{i}.prob_uc(isInTimeWin{iTimeWin}) = pu;
        tempStats{i}.h_uc(isInTimeWin{iTimeWin}) = hu;
    end
end

genTimeStr = 'tr';
for i = 1:numel(anNames)
    pName = sprintf('p_%s_%s',genTimeStr,anNames{i});
    hName = sprintf('h_%s_%s',genTimeStr,anNames{i});
    stName = sprintf('st_%s_%s',genTimeStr,anNames{i});
    stats.(pName) = tempStats{i}.prob;
    stats.(hName) = tempStats{i}.h;
    stats.(stName) = tempStats{i}.stat;
    stats.(sprintf('pu_%s_%s',genTimeStr,anNames{i})) = tempStats{i}.prob_uc;
    stats.(sprintf('hu_%s_%s',genTimeStr,anNames{i})) = tempStats{i}.h_uc;
end

end

% for i = 1:numel(anNames)
%     % Choosing variables for the actual analysis
%     [~,actVarIdx] = ismember(anVarNames{i},dataFieldNames);
%     temp = indivData(actVarIdx,:);
%     % In case of a one sample test generate matching zero data
%     if size(temp,1) == 1
%         temp(2,:) = cellfun(@(x) zeros(size(x)),temp,'UniformOutput',false);
%     end
%     temp = temp(:)';
%     data = cell(size(temp));
%     [tempStats{i}.prob,tempStats{i}.stat] = deal(NaN(size(time)));
%     % I'm using fieldtrip's high level statistical functions to perform the
%     % cluster based permutation testing. Therefore the data must be 
%     % converted to fieldtrip timelock datasturcture. See the 
%     % documentation of ft_datatype_timelock for details.
    
%     for iTimeWin = 1:numel(isInTimeWin)
%         for j = 1:size(temp,2)
%             data{j}.dimord = 'chan_time';
%             data{j}.avg = temp{j}(isInTimeWin{iTimeWin});
%             if iscolumn(data{j}.avg), data{j}.avg = data{j}.avg'; end
%             data{j}.label = {'foo'};
%             data{j}.time = time(isInTimeWin{iTimeWin});
%         end
%         % Creating the design
%         ivars = repmat(ff,1,nSubj);
%         uvar = repmat(1:nSubj,nCond,1);
%         uvar = uvar(:)';
%         design = cat(1,uvar,ivars);
%         % Configuration structure
%         % settings for ft_timelockstatistics
%         cfg.method           = 'montecarlo';
%         cfg.avgoverchan      = 'yes';
%         cfg.parameter        = 'avg';
%         % settings for ft_statistics_montecarlo
%         cfg.design           = design;
%         cfg.numrandomization = 1000;
%         cfg.correctm         = 'cluster';
%         cfg.alpha            = 0.05;
%         cfg.tail             = 1;
%         cfg.ivar             = 2;
%         cfg.uvar             = 1;
%         cfg.randomseed       = 'yes';
%         cfg.feedback         = 'text';
%         cfg.clusterstatistic = 'maxsum';
%         cfg.clusterthreshold = 'parametric';
%         cfg.clusteralpha     = 0.05;
%         cfg.clustertail      = 1;
%         cfg.statistic        = 'ft_statfun_depsamplesT';
%         % Computing statistics 
%         rawStats = ft_timelockstatistics(cfg,data{:});
%         tempStats{i}.prob(isInTimeWin{iTimeWin}) = rawStats.prob';
%         tempStats{i}.stat(isInTimeWin{iTimeWin}) = rawStats.stat';
%     end
% end
