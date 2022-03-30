function logBF_table = figure_PCM_EEG_prepost(dataPath,varargin)

% Parsing input
p = inputParser;

validThirdModels = {'Pmf','Rt'};

addRequired(p,'dataPath',@ischar);
addParameter(p,'thirdModel','',@(x) ismember(x,validThirdModels));
addParameter(p,'pltYvalTxt',false,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'pltIndiv',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'connectSubj',false,@(x) validateattributes(x,{'logical'},{'scalar'}));

parse(p,dataPath,varargin{:});

dataPath = p.Results.dataPath;
thirdModel = p.Results.thirdModel;
pltYvalTxt = p.Results.pltYvalTxt;
pltIndiv = p.Results.pltIndiv;
connectSubj = p.Results.connectSubj;

if isempty(thirdModel)
    dataRand = load(dataPath);
    out = dataRand.out;
    M = dataRand.M;
    nTimeWin = size(out,1);
    plotOrder = repmat({2:5},1,nTimeWin);
    modelNames = {{'Sp+Dec','SpRecal+Dec','Sp+DecRecal','SpRecal+DecRecal'}};
    if pltIndiv
        yBreaks = {[-9.5,-2.5],[],[],[]};
        yLims = {[-9.7,1],[],[],[]};
    else
        yBreaks = {[-9.5,-0.5],[],[],[]};
        yLims = {[-9.7,0.4],[],[],[]};
    end
elseif strcmp(thirdModel,'Pmf')
    dataRand = load(dataPath);
    out = dataRand.out;
    M = dataRand.M;
    nTimeWin = size(out,1);
    plotOrder = repmat({2:9},1,nTimeWin);
    modelNames = {{'Sp+Dec+Pmf','SpRecal+Dec+Pmf','Sp+DecRecal+Pmf',...
                   'Sp+Dec+PmfRecal','Sp+DecRecal+PmfRecal',...
                   'SpRecal+Dec+PmfRecal','SpRecal+DecRecal+Pmf',...
                   'SpRecal+DecRecal+PmfRecal'}};
    yBreaks = {[],[],[],[]};
    yLims = {[],[],[],[]};
end

nPlot = numel(plotOrder);
twNames = repmat({'50-150 ms','150-250 ms','250-350 ms','350-450 ms'},1,2);
colors = [.7 0 0;...      % red
          0 0 .7;...      % blue
          .9 .6 0;...     % orange
          0 0.6 0.6;...   % cyan
          0.5 0 0.5];     % purple

figure();
set(gcf, 'Position', [0 0 1350 250], 'PaperPositionMode', 'auto');
[mean_logBF,sem_logBF,indiv_logBF] = deal(cell(1,size(modelNames{1},2)));
vfcn = @(x) std(x)/sqrt(size(x,1));
for iPlot = 1:nPlot
    subplot(1,nPlot,iPlot);
    Tgroup = out{iPlot,1};
    Tcross = out{iPlot,2};
    nullModel = plotOrder{iPlot}(1);
    T = pcm_plotModelLikelihood(Tcross,M,'Nnull',nullModel,...
        'mindx',plotOrder{iPlot},...
        'style','bar'); %,'colors',repmat({[1,1,1]},1,nModel)
    % Plot mean values above bars
    if pltYvalTxt
        Y = mean(T.likelihood_norm(:,plotOrder{iPlot}))';
        ycoord = Y;
        ycoord(ycoord < 0) = 0;
        text(1:length(Y),ycoord,num2str(Y,'%0.1f'),'HorizontalAlignment','center',...
            'VerticalAlignment','bottom','FontSize',7,'Margin',0.01,...
            'BackgroundColor','w');
    end
    if pltIndiv
        data_indiv = T.likelihood_norm(:,plotOrder{iPlot});
        hold on;
        x = repmat((1:size(data_indiv,2))-0.25,size(data_indiv,1),1);
        if connectSubj
            plot(x',data_indiv','-o','MarkerSize',5);
        else
            scatter(x(:),data_indiv(:),20,'MarkerEdgeColor','k');
        end
    end
    set(gca,'XTickLabels',[]);
    title(twNames{iPlot});
    if iPlot == 1
        ylabel('logBF');
    else
        ylabel('');
    end
    if iPlot == nPlot
        h = findobj(gca,'Type','Bar');
        legend(flipud(h),modelNames{1});
    end
    
    if ~isempty(yLims{iPlot}), ylim(yLims{iPlot}); end
    
    if ~isempty(yBreaks{iPlot})
        breakyaxis(yBreaks{iPlot},0.005);
    else
        
    end
    
    % Saving data for table
    modelSel = [plotOrder{iPlot},size(T.likelihood_norm,2)];
    mean_logBF{iPlot} = mean(T.likelihood_norm(:,modelSel));
    sem_logBF{iPlot} = vfcn(T.likelihood_norm(:,modelSel));
    indiv_logBF{iPlot} = T.likelihood_norm(:,modelSel);
end
% Arrange BFs into a table
modelNames = cat(2,modelNames{1},{'freechol'})';
modelIdx = repmat(modelNames,nTimeWin,1);
twIdx = repmat(twNames(1:nTimeWin),numel(modelNames),1);
twIdx = twIdx(:);
tempMean = cellfun(@transpose,mean_logBF,'UniformOutput',false);
tempMean = cat(1,tempMean{:});
tempSem = cellfun(@transpose,sem_logBF,'UniformOutput',false);
tempSem = cat(1,tempSem{:});
str = arrayfun(@(f1,f2) sprintf('%.1f +/- %.1f',f1,f2),tempMean,tempSem,...
               'UniformOutput',false);
tempIndiv = cellfun(@transpose,indiv_logBF,'UniformOutput',false);
tempIndiv = cat(1,tempIndiv{:});
logBF_table = table(twIdx,modelIdx,tempMean,tempSem,str,tempIndiv,'VariableNames',...
                    {'timewin','model','mean','sem','str','indiv'});

end