function hFig = showRecalIndex(obj,varargin)
% Method for plotting estimated AV model weigths
% 
% USAGE: 
%   hFig = showRecalIndex(obj)
%   hFig = showRecalIndex(obj,'Name',Value)
% INPUT:
%   Required:
%       obj (object): mvpares object
%   'Name'-Value arguments:
%       addStats (logical): wheter to add markers indicating significant
%           effects. Default: true.
%       genTime (string): indicating the generalization time. Possible 
%           values: 'tr','tr_x_tr' for trainin time and training time by 
%           training time respectively. Default: 'tr'.
%       inference (string): group level inference 'rfx' for random
%           effects analysis, 'ffx' for fixed effects analysis
%       kind (string): deg or bin
%       plotError (logical): wheter to add error bars (takes effect only if
%           time courses are plotted (e.g. genTime is 'tr'). Default: true
%       smooth (logical): whether to smooth the time course of AV estimates
%       cust (struct): various settings to customize the plot's
%           appearence
%           Possible fields: supTitle, degLims
% OUTPUT:
%   hFig (scalar/vector): figure handles for plotted figures  

% Copyright(C) 2017, Mate Aller
% allermat@gmail.com

% Parsing input
p = inputParser;
validGenTimes = {'tr','tr_x_tr'};
validInference = {'rfx','ffx'};
validKinds = {'deg','bin'};
addRequired(p,'obj');
addParameter(p,'addStats',true,@(x)validateattributes(x,{'logical'},...
    {'scalar'}));
addParameter(p,'genTime','tr',@(x) any(validatestring(x,validGenTimes)));
addParameter(p,'inference','rfx',@(x) any(validatestring(x, ...
                                                  validInference)));
addParameter(p,'kind','deg',@(x) any(validatestring(x,validKinds)));
addParameter(p,'plotError',true,@(x)validateattributes(x,{'logical'},...
    {'scalar'}));
addParameter(p,'smooth',false,@(x)validateattributes(x,{'logical'},...
    {'scalar'}));
addParameter(p,'cust',[],@isstruct);

parse(p,obj,varargin{:});

obj = p.Results.obj;
addStats = p.Results.addStats;
genTime = p.Results.genTime;
inference = p.Results.inference;
kind = p.Results.kind;
plotError = p.Results.plotError;
smooth = p.Results.smooth;
cust = p.Results.cust;

RI = obj.getRecalIndex('genTime',genTime,'inference',inference,'smooth',smooth);
if isempty(RI)
    warning('mvpares:showRecalIndex:reqestedDataNotPresent',...
        ['Couldn''t find recalibration index in the dataset, ' ...
         'returning.']);
    hFig = [];
    return;
end
info = obj.getInfo;

if strcmp(obj.level,'subj')
    titleExtStr = [', subject ',info.subID];
else
    titleExtStr = [', group level, ',inference];
end

if strcmp(kind,'deg')
    idxName = 'RI';
else
    idxName = 'RIbin';
end

if addStats
    stats = obj.getStats('recalIndex','genTime',genTime,'smooth',smooth);
    if isempty(stats)
        warning('mvpares:showRecalIndex:reqestedDataNotPresent',...
            ['Couldn''t find statistics for recalibration index in the dataset, ',...
            'plotting without statistics.']);
        addStats = false;
    else
        statFields = fieldnames(stats);
        stats = struct2cell(stats);
    end
else
    stats = [];
end

time = obj.getTrTimePoints;
sTime = obj.plotTimeWin(1);
eTime = obj.plotTimeWin(2);
if isfield(cust,'degLims') && ~isempty(cust.degLims)
    degLims = cust.degLims;
else
    degLims = [-1.5,2.5];
end

hFig = figure();
det.hFig = hFig;
if strcmp(kind,'deg')
    det.dataUnit = 'RI (degree)';
else
    det.dataUnit = 'RI (%)';
end
switch genTime
  case 'tr'
    data = RI.(idxName);
    if strcmp(obj.level,'subj')
        error = zeros(size(data));
    else
        error = RI.([idxName,'_err']);
    end
    % Adding statistics if appropriate
    if addStats
        idx = ~cellfun(@isempty,regexp(statFields,...
                                       ['h_',idxName],'once'));
        actH = stats{idx};
        clusterIdx = mvpares.findClusters(actH);
        if ~isempty(clusterIdx)
            det.addRect = arrayfun(@(x) time(x),clusterIdx);
        else
            det.addRect = [];
        end
    end
    det.lineProp = {'Color','k','LineWidth',1.5,'LineStyle','-'};
    det.plotError = plotError;
    det.title = ['Recalibration index',titleExtStr];
    % det.yLim = degLims;
    det.xLim = [sTime,eTime];
    hFig = mvpares.plotTimeSeriesData(time,data,error,det);
    % set(hFig,'Units','normalized','Position',[1/4,1/4,1/2,1/2],...
                % 'PaperPositionMode','auto');
  case 'tr_x_tr'
    hFig = [];
    % condsOfInterest = {'R_D_a','R_D_v','R_d_a','R_d_v',...
    %                    'r_D_a','r_D_v','r_d_a','r_d_v'};
    % titles = {'VR+, D+, A','VR+, D+, V','VR+, D-, A','VR+, D-, V',...
    %           'VR-, D+, A','VR-, D+, V','VR-, D-, A','VR-, D-, V'};
    % nColFig = 4;
    % nRowFig = ceil(numel(condsOfInterest)/nColFig);
    % for i = 1:numel(condsOfInterest)
    %     data = weights.(strcat(condsOfInterest{i},'_wav'));
    %     det.hAxes = subplot(nRowFig,nColFig,i);
    %     det.title = titles{i};
    %     det.cLim = degLims;
    %     det.xLim = [sTime,eTime];
    %     mvpares.plotTimeByTimeMatrix(time,data,det);
    % end
    % set(hFig,'Units','normalized','Position',[1/8,1/4,3/4,0.6],...
    %          'PaperPositionMode','auto');
    % temp = findobj(hFig,'Type','axes');
    % hAxes = temp(1);
    % posOutAxes = get(hAxes,'OuterPosition');
    % posOutCbar = [posOutAxes(1)+(posOutAxes(3)),posOutAxes(2),...
    %               (posOutAxes(3)*0.3),posOutAxes(4)];
    % hCbar = colorbar('OuterPosition',posOutCbar);
    % ylabel(hCbar,'W_A_V (degree)');
    % if isfield(cust,'supTitle')
    %     suplabel(cust.supTitle,'t',[0.08,0.08,0.84,0.88]);
    % else
    %     suplabel('W_A_V 3-way interaction','t',[0.08,0.08,0.84,0.88]);
    % end
end


end