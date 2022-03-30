function logBF_table = figure_PCM_fMRI_prepost(dataPath,varargin)

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

% Load data
if isempty(thirdModel)
    % Data from the random effects analysis
    dataRand = load(dataPath);
    if pltIndiv
        yBreaks = {[-9.3,-1.5],[12,35],[],[],[]};
        yLims = {[-9.5,1.6],[-1,36],[-1,22],[-5,60],[-2,30]};
    else
        yBreaks = {[-7.9,-0.25],[9,67],[16,30],[42,115],[21,45]};
        yLims = {[-8.1,1.5],[0,68],[0,31],[0,117],[0,47]};
    end
    nROI = size(dataRand.out,1);
    plotOrder = repmat({2:5},1,nROI);
    modelNames = {{'Sp+Dec','SpRecal+Dec','Sp+DecRecal','SpRecal+DecRecal'}};
elseif strcmp(thirdModel,'Pmf')
    % Data from the random effects analysis
    dataRand = load(dataPath);
    yBreaks = {[],[],[],[],[],[],[],[]};
    yLims = {[],[],[],[],[],[],[],[]};
    nROI = size(dataRand.out,1);
    plotOrder = repmat({2:9},1,nROI);
    modelNames = {{'Sp+Dec+Pmf','SpRecal+Dec+Pmf','Sp+DecRecal+Pmf',...
                   'Sp+Dec+PmfRecal','Sp+DecRecal+PmfRecal',...
                   'SpRecal+Dec+PmfRecal','SpRecal+DecRecal+Pmf',...
                   'SpRecal+DecRecal+PmfRecal'}};
elseif strcmp(thirdModel,'Rt')
    error('Not yet implemented!');
end

roiNames = repmat({'HG','hA','IPL','IPS','FEF'},1,2);
[mean_logBF,sem_logBF,indiv_logBF] = deal(cell(1,numel(plotOrder)));
vfcn = @(x) std(x)/sqrt(size(x,1));
iPlot = 1;
for iFig = 1:1
    figure();
    set(gcf, 'Position', [0 0 1350 250], 'PaperPositionMode', 'auto');
    for iROI = 1:nROI
        subplot(1,nROI,iROI);
        out = dataRand.out;
        M = dataRand.M;
        nullModel = plotOrder{iFig,iROI}(1);
        Tgroup = out{iROI,1};
        Tcross = out{iROI,2};
        T = pcm_plotModelLikelihood(Tcross,M,'Nnull',nullModel,...
            'mindx',plotOrder{iFig,iROI},...
            'style','bar'); %,'colors',repmat({[1,1,1]},1,nModel)
        % Plot mean values above bars
        if pltYvalTxt
            Y = mean(T.likelihood_norm(:,plotOrder{iFig,iROI}))';
            ycoord = Y;
            ycoord(ycoord < 0) = 0;
            text(1:length(Y),ycoord,num2str(Y,'%0.1f'),'HorizontalAlignment','center',...
                'VerticalAlignment','bottom','FontSize',7,'Margin',0.01,...
                'BackgroundColor','w');
        end
        if pltIndiv
            data_indiv = T.likelihood_norm(:,plotOrder{iFig,iROI});
            hold on;
            x = repmat((1:size(data_indiv,2))-0.25,size(data_indiv,1),1);
            if connectSubj
                plot(x',data_indiv','-o','MarkerSize',5);
            else
                scatter(x(:),data_indiv(:),20,'MarkerEdgeColor','k');
            end
        end
        set(gca,'XTickLabels',[]);
        title(roiNames{iROI});
        if iROI == 1
            ylabel('logBF');
        else
            ylabel('');
        end
        if iROI == nROI
            h = findobj(gca,'Type','Bar');
            legend(flipud(h),modelNames{iFig});
        end
        if ~isempty(yLims{iFig,iROI}), ylim(yLims{iFig,iROI}); end
        if ~isempty(yBreaks{iFig,iROI})
            breakyaxis(yBreaks{iFig,iROI},0.005);
        else
            
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
modelNames = cat(2,modelNames{1},{'freechol'})';
modelIdx = repmat(modelNames,nROI,1);
roiIdx = repmat(roiNames(1:nROI),numel(modelNames),1);
roiIdx = roiIdx(:);
tempMean = cellfun(@transpose,mean_logBF,'UniformOutput',false);
tempMean = cat(1,tempMean{:});
tempSem = cellfun(@transpose,sem_logBF,'UniformOutput',false);
tempSem = cat(1,tempSem{:});
% tempMean = arrayfun(@(x) cat(2,mean_logBF{1,x},mean_logBF{2,x})',1:nROI,'UniformOutput',false);
% tempMean = cat(1,tempMean{:});
% tempSem = arrayfun(@(x) cat(2,sem_logBF{1,x},sem_logBF{2,x})',1:nROI,'UniformOutput',false);
% tempSem = cat(1,tempSem{:});
str = arrayfun(@(f1,f2) sprintf('%.1f +/- %.1f',f1,f2),tempMean,tempSem,...
               'UniformOutput',false);
tempIndiv = cellfun(@transpose,indiv_logBF,'UniformOutput',false);
tempIndiv = cat(1,tempIndiv{:});
logBF_table = table(roiIdx,modelIdx,tempMean,tempSem,str,tempIndiv,...
                    'VariableNames',{'roi','model','mean','sem','str','indiv'});
end