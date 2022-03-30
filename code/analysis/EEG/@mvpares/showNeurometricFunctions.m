function hFig = showNeurometricFunctions(obj,varargin)
% Method for plotting neurometric functions
%
% USAGE:
%   hFig = showNeurometricFunctions(obj)
%   hFig = showNeurometricFunctions(obj,'Name',Value)
% INPUT:
%   Required:
%       obj (object): mvpares object
%   'Name'-Value arguments:
%       inference (string): group level inference 'rfx' for random
%           effects analysis, 'ffx' for fixed effects analysis
%       timePoint (vector):
%       in (structure):
%
% OUTPUT:
%   hFig (scalar): figure handle for plotted figure

% Copyright(C) 2017, Mate Aller
% allermat@gmail.com

% Parsing input, checking matlab
p = inputParser;

validGenTimes = {'tr','tr_x_tr'};
validInference = {'rfx','ffx'};
validTreatments = {'individual','collapse'};
addRequired(p,'obj');
addParameter(p,'inference','rfx',@(x) any(validatestring(x,validInference)));
addParameter(p,'timePoint',[],@(x) validateattributes(x,{'numeric'},{'finite'}));
addParameter(p,'treatTimePoints','individual',@(x) any(validatestring(x,validTreatments)));
addParameter(p,'in',[],@isstruct);

parse(p,obj,varargin{:});

obj = p.Results.obj;
inference = p.Results.inference;
timePoint = p.Results.timePoint;
treatTimePoints = p.Results.treatTimePoints;
in = p.Results.in;

% Loading necessary data
if strcmp(obj.level,'subj')
    nmf = obj.getNeuroMetricFunctions('genTime','tr');
else
    nmf = obj.getNeuroMetricFunctions('genTime','tr','inference',inference);
end

trTimePoints = obj.getTrTimePoints;
if isempty(nmf)
    hFig = [];
    warning('mvpares:showNeurometricFunctions:reqestedDataNotPresent',...
            'Couldn''t find neurometric functions in the dataset, returning.');
    return;
end

if ~isempty(timePoint)
    timePointIdx = obj.indTrTimePoint(timePoint);
    if isempty(timePointIdx)
        hFig = [];
        warning('mvpares:showNeurometricFunctions:invalidInput',...
                ['The requested time point is not present in the dataset', ...
                 ', returning. ']);
    end
else
    timePointIdx = 1:size(trTimePoints,1);
end

% In case the specified time window is to be treated as collapsed
% we need to fit the neurometric functions
if strcmp(treatTimePoints,'collapse')
    PF = nmf.PF;
    sl = nmf.stimLevels;
    tempNmf = nmf;
    fNames = fieldnames(nmf);
    fieldsToRm = fNames(cellfun(@isempty,regexp(fNames,'_pctr$')));
    tempNmf = rmfield(tempNmf,fieldsToRm);
    tempNmf = structfun(@(x) mean(x(:,timePointIdx,:),2),tempNmf,'UniformOutput',false);
    tempPct = structfun(@(x) squeeze(mat2cell(x,size(x,1),size(x,2),ones(size(x,3),1))),...
                        tempNmf,'UniformOutput',false);
    % Fitting parameters
    opt.gamma = 0.2;
    opt.lambda = 0.2;
    fun_fit = @(x) fit_PF2(PF,sl,(x*100),ones(size(sl))*100,opt);
    tempPct = structfun(@(x) cellfun(fun_fit,x,'UniformOutput',false),...
                        tempPct,'UniformOutput',false);
    tempPct = structfun(@(x) cell2mat(shiftdim(cellfun(@transpose,x,'UniformOutput',false),-2)),...
                        tempPct,'UniformOutput',false);
    fNames = fieldnames(tempPct);
    newFnames = strrep(fNames,'_pctr','_pv');
    for i = 1:numel(newFnames)
        [tempNmf.(newFnames{i})] = tempPct.(fNames{i});
    end
    tempNmf.PF = PF;
    tempNmf.stimLevels = sl;
    nmf = tempNmf;
    % Updating timePointIdx as it is no longer a vector
    timePointIdx = 1;
end

nmGroupingLevels = {'pre','post_L','post_R'};
lineProp =  cat(1,...
                {'Color',[0,0.6,0],'LineWidth',2},...
                {'Color',[0,0,0.6],'LineWidth',2},...
                {'Color',[0.6,0,0],'LineWidth',2});
markerProp = cat(1,...
                 {'ko','MarkerSize',10,'MarkerFaceColor',[0,0.6,0]},...
                 {'ko','MarkerSize',10,'MarkerFaceColor',[0,0,0.6]},...
                 {'ko','MarkerSize',10,'MarkerFaceColor',[0.6,0,0]});
for iTimePoint = 1:numel(timePointIdx)
    hCurve = [];
    legends = {'pre-test','L adaptation','R adaptation'};
    hFig(iTimePoint) = figure(); %#ok
    parentAxes = axes('parent',hFig(iTimePoint));
    hold on;
    for i = 1:numel(nmGroupingLevels)
        if ~isfield(nmf,[nmGroupingLevels{i},'_pv'])
            continue;
        end
        actPv = nmf.([nmGroupingLevels{i},'_pv']);
        actPctr = nmf.([nmGroupingLevels{i},'_pctr']);        

        if strcmp(obj.level,'group') && strcmp(inference,'rfx')
            
            subjDim = numel(size(actPctr));
            actPctr = mean(actPctr,subjDim);
            % Plot data
            plot(parentAxes,nmf.stimLevels,actPctr(:,timePointIdx(iTimePoint)),markerProp{i,:});

            % % Plot PF fit
            StimLevelsFineGrain = min(nmf.stimLevels):max(nmf.stimLevels)/1000:max(nmf.stimLevels);
            for iSubj = 1:numel(size(actPv))
            ProportionCorrectModel(:,iSubj) = nmf.PF(actPv(:, ... 
                                                  timePointIdx(iTimePoint),iSubj),StimLevelsFineGrain); %#ok
            end
            ProportionCorrectModel_mean = mean(ProportionCorrectModel,2);
            ProportionCorrectModel_sem = std(ProportionCorrectModel,0,2)./size(actPv,subjDim);
            h = shadedErrorBar(StimLevelsFineGrain,ProportionCorrectModel_mean,...
                                  ProportionCorrectModel_sem,lineProp(i,:),1); 
            hCurve(i) = h.mainLine; %#ok
            % hCurve(i) = plot(parentAxes,,); %#ok
        else
            
            % Plot data
            plot(parentAxes,nmf.stimLevels,actPctr(:,timePointIdx(iTimePoint)),markerProp{i,:});

            % % Plot PF fit
            StimLevelsFineGrain = min(nmf.stimLevels):max(nmf.stimLevels)/1000:max(nmf.stimLevels);
            ProportionCorrectModel = nmf.PF(actPv(:,timePointIdx(iTimePoint)),StimLevelsFineGrain);
            hCurve(i) = plot(parentAxes,StimLevelsFineGrain, ...
                             ProportionCorrectModel,lineProp{i,:}); %#ok
        end
    end
    % Default format figure
    set(parentAxes, 'FontSize', 14, 'Xtick', nmf.stimLevels);
    xlabel(parentAxes, 'Sound location (azimuth)', 'fontsize', 18)
    ylabel(parentAxes, '% decoded right', 'fontsize', 18)
    xlim(parentAxes, [-15,15]);
    ylim(parentAxes, [0,1])
    legends = legends(hCurve ~= 0);
    hCurve = hCurve(hCurve ~= 0);
    legend(hCurve,legends{:},'location','NorthWest');
    title(parentAxes,sprintf('%d ms',round(trTimePoints(timePointIdx(iTimePoint))*1000)));
end

end