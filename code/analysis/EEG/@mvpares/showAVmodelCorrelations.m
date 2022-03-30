function hFig = showAVmodelCorrelations(obj,var,varargin)
% Method for plotting AV model correlations
%
% USAGE:
%   hFig = showAVmodelCorrelations(obj)
%   hFig = showAVmodelCorrelations(obj,'Name',Value)
% INPUT:
%   Required:
%       obj (object): mvpares object
%       var (string): the variable whos AV weights were correlated with 
%           the AV weights of the mvpares object. Possible values:
%           'acrossTime', 'behav', 'fmri'
%   'Name'-Value arguments:
%       genTime (string): indicating the generalization time. Possible 
%           values: 'tr','tr_x_tr' for trainin time and training time by 
%           training time respectively. Default: 'tr'. 
%       smooth (logical): whether to smooth the time courses 
% OUTPUT:
%   hFig (scalar/vector): figure handles for plotted figures 
%
% Copyright(C) 2016, Mate Aller

% Parsing input
p = inputParser;
validVars = {'acrossTime','behav','fmri'};
validGenTimes = {'tr','tr_x_tr'};
addRequired(p,'obj');
addRequired(p,'var',@(x) any(validatestring(x,validVars)));
addParameter(p,'genTime','tr',@(x) any(validatestring(x,validGenTimes)));
addParameter(p,'smooth',false,@(x) validateattributes(x,{'logical'},...
    {'scalar'}));
parse(p,obj,var,varargin{:});
obj = p.Results.obj;
genTime = p.Results.genTime;
var = p.Results.var;
smooth = p.Results.smooth;

corrCoeffs = obj.getAVmodelCorrelations(var,'genTime',genTime,'smooth',smooth);

if isempty(corrCoeffs)
    hFig = [];
    warning('mvpares:showAVmodelEstimates:reqestedDataNotPresent',...
        'Couldn''t find AV model correlations in the dataset, returning.');
    return;
end

hFig = [];
time = obj.getTrTimePoints;
tWin = obj.plotTimeWin;

switch var
    case 'acrossTime'
        data = corrCoeffs.cc_acrossTime;
        error = corrCoeffs.ci_acrossTime;
        det.title = 'Across time correlation of neural AV weights';
        det.dataUnit = 'R';
        det.xLim = tWin;
        switch genTime
            case 'tr'
                det.lineProp = {'Color','k','LineWidth',1.5,'LineStyle','-'};
                hFig = mvpares.plotTimeSeriesData(time,data,error,det);
            case 'tr_x_tr'
                det.cLim = [-1,1];
                det.addCbar = true;
                det.xLab = 'time (ms)';
                det.yLab = 'time (ms)';
                hFig = mvpares.plotTimeByTimeMatrix(time,data,det);
        end
    case 'behav'
        data = corrCoeffs.cc_behav;
        error = corrCoeffs.ci_behav;
        det.title = 'Correlation between neural and behavioural AV weights';
        det.dataUnit = 'R';
        det.xLim = tWin;
        switch genTime
            case 'tr'
                det.lineProp = {'Color','k','LineWidth',1.5,'LineStyle','-'};
                det.yLim = [min(data),1];
                hFig = mvpares.plotTimeSeriesData(time,data,error,det);
            case 'tr_x_tr'
                det.cLim = [-1,1];
                det.addCbar = true;
                hFig = mvpares.plotTimeByTimeMatrix(time,data,det);
        end
    case 'fmri'
        fields = fieldnames(corrCoeffs);
        fmriCorrNames = fields(~cellfun(@isempty,regexp(fields,'^cc_fmri','once')));
        % Reordering according to the hierarchy
        fmriCorrNames = fmriCorrNames([4:7,2:3,1,8]);
        temp = regexp(fmriCorrNames,'^cc_fmri_(.*)','tokens');
        rois = [temp{:}]';
        hFig = figure();
        det.hFig = hFig;
        nColFig = 4;
        nRowFig = ceil(numel(fmriCorrNames)/nColFig);
        switch genTime
            case 'tr'
                for i = 1:numel(fmriCorrNames)
                    data = corrCoeffs.(fmriCorrNames{i});
                    error = corrCoeffs.(strrep(fmriCorrNames{i},'cc_','ci_'));
                    det.hAxes = subplot(nRowFig,nColFig,i);
                    det.title = strrep(rois{i},'_','\_');
                    det.dataUnit = 'R';
                    det.lineProp = {'Color','k','LineWidth',1.5,'LineStyle','-'};
                    det.yLim = [-0.2,1];
                    det.xLim = tWin;
                    mvpares.plotTimeSeriesData(time,data,error,det);
                end
                set(hFig,'Units','normalized','Position',[1/8,1/4,3/4,0.6],...
                    'PaperPositionMode','auto');
                suplabel('Correlation between EEG and fMRI AV weights','t');
            case 'tr_x_tr'
                for i = 1:numel(fmriCorrNames)
                    data = corrCoeffs.(fmriCorrNames{i});
                    det.hAxes = subplot(nRowFig,nColFig,i);
                    det.title = strrep(rois{i},'_','\_');
                    det.cLim = [-1,1];
                    det.xLim = tWin;
                    mvpares.plotTimeByTimeMatrix(time,data,det);
                end
                set(hFig,'Units','normalized','Position',[1/8,1/4,3/4,0.6],...
                    'PaperPositionMode','auto');
                temp = findobj(hFig,'Type','axes');
                hAxes = temp(1);
                posOutAxes = get(hAxes,'OuterPosition');
                posOutCbar = [posOutAxes(1)+(posOutAxes(3)),posOutAxes(2),...
                    (posOutAxes(3)*0.3),posOutAxes(4)];
                hCbar = colorbar('OuterPosition',posOutCbar);
                ylabel(hCbar,'R');
                suplabel('Correlation between EEG and fMRI AV weights','t');
        end
end

end