function logBF_table = figure_PCM_fMRI_to_EEG_prepost(dataPath,varargin)
% Script for plotting figure with the winning models from the original
% model fitting and see if the PMF modal adds any evidence

% Parsing input
p = inputParser;

addRequired(p,'dataPath',@ischar);
addParameter(p,'pltYvalTxt',false,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'pltIndiv',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'connectSubj',false,@(x) validateattributes(x,{'logical'},{'scalar'}));

parse(p,dataPath,varargin{:});

dataPath = p.Results.dataPath;
pltYvalTxt = p.Results.pltYvalTxt;
pltIndiv = p.Results.pltIndiv;
connectSubj = p.Results.connectSubj;

load(dataPath,'out','M');

nPlot = size(out,1);
twNames = {'50-150 ms','150-250 ms','250-350 ms','350-450 ms'};
nTw = numel(twNames);
if pltIndiv
    yBreaks = {[24,47],[75,120],[],[]};
    yLims = {[],[-5,130],[0,160],[-15,95]};
else
    yBreaks = {[20,47],[50,120],[100,145],[60,90]};
    yLims = {[],[0,130],[0,160],[0,95]};
end
xTickLabels = {'HG','hA','IPL','IPS','FEF'};
subjCols = [.7 0 0;...      % red
           0 0 .7;...      % blue
           .9 .6 0;...     % orange
           0 0.6 0.6;...   % cyan
           0.5 0 0.5];     % purple
modelsToPlot = 2:6;
[mean_logBF,sem_logBF,indiv_logBF] = deal(cell(1,nPlot));
vfcn = @(x) std(x)/sqrt(size(x,1));
figure();
set(gcf, 'Position', [0 0 1350 250], 'PaperPositionMode', 'auto');
for iPlot = 1:nPlot
    subplot(1,nPlot,iPlot);
    Tgroup = out{iPlot,1};
    Tcross = out{iPlot,2};
    nullModel = 2;
    T = pcm_plotModelLikelihood(Tcross,M,'Nnull',nullModel,...
        'mindx',modelsToPlot,...
        'style','bar'); %,'colors',repmat({[1,1,1]},1,nModel)
    % Plot mean values above bars
    if pltYvalTxt
        Y = mean(T.likelihood_norm(:,modelsToPlot))';
        ycoord = Y;
        ycoord(ycoord < 0) = 0;
        text(1:length(Y),ycoord,num2str(Y,'%0.1f'),'HorizontalAlignment','center',...
            'VerticalAlignment','bottom','FontSize',7,'Margin',0.01,...
            'BackgroundColor','w');
    end
    if pltIndiv
        data_indiv = T.likelihood_norm(:,modelsToPlot);
        hold on;
        x = repmat((1:size(data_indiv,2))-0.25,size(data_indiv,1),1);
        if connectSubj
            plot(x',data_indiv','-o','MarkerSize',5);
        else
            scatter(x(:),data_indiv(:),20,'MarkerEdgeColor','k');
        end
    end
    set(gca,'XTickLabels',[]);
    %         set(gca,'XTickLabels',xTickLabels{iPlot});
    %         xtickangle(45);
    title(twNames{iPlot});
    if iPlot == 1
        ylabel('logBF');
    else
        ylabel('');
    end
    if iPlot == nPlot
        h = findobj(gca,'Type','Bar');
        legend(flipud(h),xTickLabels);
    end
    
    if ~isempty(yLims{iPlot}), ylim(yLims{iPlot}); end
    
    if ~isempty(yBreaks{iPlot})
        breakyaxis(yBreaks{iPlot},0.005);
    else
        
    end
    
    % Saving data for table
    modelSel = [modelsToPlot,size(T.likelihood_norm,2)];
    mean_logBF{iPlot} = mean(T.likelihood_norm(:,modelSel))';
    sem_logBF{iPlot} = vfcn(T.likelihood_norm(:,modelSel))';
    indiv_logBF{iPlot} = T.likelihood_norm(:,modelSel)';
end
% Arrange BFs into a table
modelNames = cat(2,xTickLabels,{'freechol'})';
modelIdx = repmat(modelNames,nTw,1);
twIdx = repmat(twNames,numel(modelNames),1);
twIdx = twIdx(:);
tempMean = cat(1,mean_logBF{:});
tempSem = cat(1,sem_logBF{:});
tempIndiv = cat(1,indiv_logBF{:});
str = arrayfun(@(f1,f2) sprintf('%.1f +/- %.1f',f1,f2),tempMean,tempSem,...
               'UniformOutput',false);
logBF_table = table(twIdx,modelIdx,tempMean,tempSem,str,tempIndiv,'VariableNames',...
                    {'timewin','model','mean','sem','str','indiv'});
end