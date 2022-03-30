function logBF_table = figure_PCM_EEG_pre_winning_pmf(dataPath,varargin)

% Parsing input
p = inputParser;

addRequired(p,'dataPath',@ischar);
addParameter(p,'pltYvalTxt',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
parse(p,dataPath,varargin{:});

dataPath = p.Results.dataPath;
pltYvalTxt = p.Results.pltYvalTxt;

dataRand = load(dataPath);
out = dataRand.out;
M = dataRand.M;
nTimeWin = size(out,1);
plotOrder = {[3,5],[3,5],[3,5],[3,5]};
modelNames = {{'Sp+Dec','Sp+Dec+Pmf'},...
              {'Sp+Dec','Sp+Dec+Pmf'},...
              {'Sp+Dec','Sp+Dec+Pmf'},...
              {'Sp+Dec','Sp+Dec+Pmf'}};
yBreaks = {[],[],[],[]};
yLims = {[],[],[],[]};
colors = {[.7 0 0],...      % red
          [0 0 .7],...      % blue
          [.9 .6 0],...     % orange
          [0 0.6 0.6],...   % cyan
          [0.5 0 0.5],...   % purple
          [0.2 0.6 0.2]};   % green
modelColors = {{colors{3},colors{6}},...
               {colors{3},colors{6}},...
               {colors{3},colors{6}},...
               {colors{3},colors{6}}};

nPlot = numel(plotOrder);
twNames = repmat({'50-150 ms','150-250 ms','250-350 ms','350-450 ms'},1,2);
figure();
set(gcf, 'Position', [0 0 1350 350], 'PaperPositionMode', 'auto');
[mean_logBF,sem_logBF,indiv_logBF] = deal(cell(1,size(modelNames{1},2)));
vfcn = @(x) std(x)/sqrt(size(x,1));
for iPlot = 1:nPlot
    subplot(1,nPlot,iPlot);
    Tgroup = out{iPlot,1};
    Tcross = out{iPlot,2};
    nullModel = 1; %plotOrder{iPlot}(1);
    T = pcm_plotModelLikelihood(Tcross,M,'Nnull',nullModel,...
        'mindx',plotOrder{iPlot},...
        'style','bar','colors',modelColors{iPlot});
    % Plot mean values above bars
    if pltYvalTxt
        Y = mean(T.likelihood_norm(:,plotOrder{iPlot}))';
        ycoord = Y;
        ycoord(ycoord < 0) = 0;
        text(1:length(Y),ycoord,num2str(Y,'%0.1f'),'HorizontalAlignment','center',...
            'VerticalAlignment','bottom','FontSize',7,'Margin',0.01,...
            'BackgroundColor','w');
    end
    set(gca,'XTick',1:numel(plotOrder{iPlot}),...
            'XTickLabel',modelNames{iPlot},'XTickLabelRotation',45);
    title(twNames{iPlot});
    if iPlot == 1
        ylabel('logBF');
    else
        ylabel('');
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
modelIdx = cellfun(@(c) cat(2,c,{'freechol'})',modelNames,'UniformOutput',false);
modelIdx = cat(1,modelIdx{:});
twIdx = repmat(twNames(1:nTimeWin),numel(modelNames{1})+1,1);
twIdx = twIdx(:);
tempMean = cellfun(@transpose,mean_logBF,'UniformOutput',false);
tempMean = cat(1,tempMean{:});
tempSem = cellfun(@transpose,sem_logBF,'UniformOutput',false);
tempSem = cat(1,tempSem{:});
str = arrayfun(@(f1,f2) sprintf('%.1f +/- %.1f',f1,f2),tempMean,tempSem,...
               'UniformOutput',false);
tempIndiv = cellfun(@transpose,indiv_logBF,'UniformOutput',false);
tempIndiv = cat(1,tempIndiv{:});
logBF_table = table(twIdx,modelIdx,tempMean,tempSem,str,tempIndiv,...
                    'VariableNames',{'timewin','model','mean','sem','str','indiv'});
end