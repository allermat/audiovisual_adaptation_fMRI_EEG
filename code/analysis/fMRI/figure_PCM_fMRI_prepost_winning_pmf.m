function logBF_table = figure_PCM_fMRI_prepost_winning_pmf(dataPath)
% Script for plotting figure with the winning models from the original
% model fitting and see if the PMF modal adds any evidence

% Parsing input
p = inputParser;

addRequired(p,'dataPath',@ischar);

parse(p,dataPath);

dataPath = p.Results.dataPath;

% Load data
% Data from the random effects analysis
dataRand = load(dataPath);

yLims = {[-0.1,1],[],[],[],[]}; % [-9.58,1]
yBreaks = {[],[],[],[],[]}; %[-9.5,-0.1]

pltYvalTxt = true;
nROI = size(dataRand.out,1);
% plotOrder = {[7:9],[2:4],[2:4],[2:4],[2:4]};
plotOrder = {[2:3],[2:3],[2:3],[2:3],[2:3]};
roiNames = {'HG','hA','IPL','IPS','FEF'};
% modelNames = {{'SpRecal+Dec','SpRecal+Dec+PmfRecal','SpRecal+Dec+Pmf'},...
modelNames = {{'SpRecal+DecRecal','SpRecal+DecRecal+PmfRecal'},...
              {'SpRecal+DecRecal','SpRecal+DecRecal+PmfRecal'},...
              {'SpRecal+DecRecal','SpRecal+DecRecal+PmfRecal'},...
              {'SpRecal+DecRecal','SpRecal+DecRecal+PmfRecal'},...
              {'SpRecal+DecRecal','SpRecal+DecRecal+PmfRecal'}};
colors = {[.7 0 0],...      % red
          [0 0 .7],...      % blue
          [.9 .6 0],...     % orange
          [0 0.6 0.6],...   % cyan
          [0.5 0 0.5],...   % purple
          [0.2 0.6 0.2]};   % green
% modelColors = {{colors{2},colors{5},colors{6}},...
modelColors = {{colors{4},colors{5},colors{6}},...
               {colors{4},colors{5},colors{6}},...
               {colors{4},colors{5},colors{6}},...
               {colors{4},colors{5},colors{6}},...
               {colors{4},colors{5},colors{6}}};
[mean_logBF,sem_logBF,indiv_logBF] = deal(cell(1,numel(plotOrder)));
vfcn = @(x) std(x)/sqrt(size(x,1));
iPlot = 1;
for iFig = 1:1
%     if iFig == 1
%         out = dataFix.out;
%         M = dataFix.M;
%         nullModel = 1;
%     else
%         out = dataRand.out;
%         M = dataRand.M;
%         nullModel = plotOrder{1,1}(1);
%     end
    out = dataRand.out;
    M = dataRand.M;
    nullModel = 5; % corresponding to spatial no recal + dec no recal
    figure();
    set(gcf, 'Position', [0 0 1350 350], 'PaperPositionMode', 'auto');
    for iROI = 1:nROI
        subplot(1,nROI,iROI);
        iROI = iROI;
        Tgroup = out{iROI,1};
        Tcross = out{iROI,2};
        T = pcm_plotModelLikelihood(Tcross,M,'Nnull',nullModel,...
                'mindx',plotOrder{iFig,iROI},...
                'style','bar','colors',modelColors{iFig,iROI});
        % Plot mean values above bars
        if pltYvalTxt
            Y = mean(T.likelihood_norm(:,plotOrder{iFig,iROI}))';
            ycoord = Y;
            ycoord(ycoord < 0) = 0;
            text(1:length(Y),ycoord,num2str(Y,'%0.1f'),'HorizontalAlignment','center',...
                'VerticalAlignment','bottom','FontSize',7,'Margin',0.01,...
                'BackgroundColor','w');
        end
        set(gca,'XTickLabels',[]);
        set(gca,'XTick',1:numel(plotOrder{iFig,iROI}),...
            'XTickLabel',modelNames{iFig,iROI},'XTickLabelRotation',45);
        title(roiNames{iROI});
        if iROI == 1
            ylabel('logBF');
        else
            ylabel('');
        end
%         if iPlot == nPlot
%             h = findobj(gca,'Type','Bar');
%             legend(flipud(h),modelNames{iFig});
%         end
        if ~isempty(yLims{iFig,iROI}), ylim(yLims{iFig,iROI}); end
        if ~isempty(yBreaks{iFig,iROI})
            breakyaxis(yBreaks{iFig,iROI},0.005);
        end
        % Saving data for table
        modelSel = [plotOrder{iFig,iROI},size(T.likelihood_norm,2)];
        mean_logBF{iPlot} = mean(T.likelihood_norm(:,modelSel));
        sem_logBF{iPlot} = vfcn(T.likelihood_norm(:,modelSel));
        indiv_logBF{iPlot} = T.likelihood_norm(:,modelSel);
        iPlot = iPlot+1;
    end
end

% Arrange BFs into a table
modelIdx = cellfun(@(c) cat(2,c,{'freechol'})',modelNames,'UniformOutput',false);
modelIdx = cat(1,modelIdx{:});
roiIdx = repmat(roiNames(1:nROI),numel(modelNames{1})+1,1);
roiIdx = roiIdx(:);
tempMean = cellfun(@transpose,mean_logBF,'UniformOutput',false);
tempMean = cat(1,tempMean{:});
tempSem = cellfun(@transpose,sem_logBF,'UniformOutput',false);
tempSem = cat(1,tempSem{:});
str = arrayfun(@(f1,f2) sprintf('%.1f +/- %.1f',f1,f2),tempMean,tempSem,...
               'UniformOutput',false);
tempIndiv = cellfun(@transpose,indiv_logBF,'UniformOutput',false);
tempIndiv = cat(1,tempIndiv{:});
logBF_table = table(roiIdx,modelIdx,tempMean,tempSem,str,tempIndiv,...
                    'VariableNames',{'roi','model','mean','sem','str','indiv'});
end