function stats = compAVmodelWeightStats(indivData,time)
% Method to compute statistics for AV weights data.

% Parsing input
p = inputParser;
addRequired(p,'indivData',@(x) validateattributes(x,{'struct'},{'vector',...
    'nonempty'}));
addRequired(p,'time',@(x) validateattributes(x,{'numeric'},{'vector',...
    'increasing'}));
parse(p,indivData,time);
indivData = p.Results.indivData;
time = p.Results.time;

if iscolumn(time), time = time'; end

% Number of subjects
nSubj = length(indivData);
% Analysis names. As for now statistics are implemented only fir main 
% effects
anNames = {'VR','Task','Disp'};
anVarNames = {{'r_wav','R_wav'},{'a_wav','v_wav'},{'d_wav','D_wav'}};
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
    temp = temp(:)';
    data = cell(size(temp));
    % I'm using fieldtrip's high level statistical functions to perform the
    % cluster based permutation testing. Therefore the data must be 
    % converted to fieldtrip timelock datasturcture. See the 
    % documentation of ft_datatype_timelock for details.
    for j = 1:size(temp,2)
        data{j}.dimord = 'chan_time';
        data{j}.avg = temp{j};
        if iscolumn(data{j}.avg), data{j}.avg = data{j}.avg'; end
        data{j}.label = {'foo'};
        data{j}.time = time;
    end
    % Creating the design
    ivars = repmat(ff,1,nSubj);
    uvar = repmat(1:nSubj,nCond,1);
    uvar = uvar(:)';
    design = cat(1,uvar,ivars);
    % Configuration structure
    % settings for ft_timelockstatistics
    cfg.method           = 'montecarlo';
    cfg.avgoverchan      = 'yes';
    cfg.parameter        = 'avg';
    % settings for ft_statistics_montecarlo
    cfg.design           = design;
    cfg.numrandomization = 1000;
    cfg.correctm         = 'cluster';
    cfg.alpha            = 0.05;
    cfg.tail             = 1;
    cfg.ivar             = 2;
    cfg.uvar             = 1;
    cfg.randomseed       = 'yes';
    cfg.feedback         = 'text';
    cfg.clusterstatistic = 'maxsum';
    cfg.clusterthreshold = 'parametric';
    cfg.clusteralpha     = 0.05;
    cfg.clustertail      = 1;
    cfg.statistic        = 'ft_statfun_circ_wwtest';
    % Computing statistics 
    tempStats{i} = ft_timelockstatistics(cfg,data{:});
end
genTimeStr = 'tr';
for i = 1:numel(anNames)
    pName = sprintf('p_%s_%s',genTimeStr,anNames{i});
    hName = sprintf('h_%s_%s',genTimeStr,anNames{i});
    stName = sprintf('st_%s_%s',genTimeStr,anNames{i});
    stats.(pName) = tempStats{i}.prob';
    stats.(hName) = stats.(pName) < 0.05;
    stats.(stName) = tempStats{i}.stat';
end

end

