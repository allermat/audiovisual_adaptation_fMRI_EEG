function pcm_plot_fitRes(fitRes,models,dataSubj,condVec,partVec,runEffect,iROI,modelsToPlot,hFig,figTitles,varargin)

% Parsing input
p = inputParser;
addRequired(p,'fitRes');
addRequired(p,'models');
addRequired(p,'dataSubj');
addRequired(p,'condVec');
addRequired(p,'partVec');
addRequired(p,'runEffect');
addRequired(p,'iROI');
addRequired(p,'modelsToPlot');
addRequired(p,'hFig');
addRequired(p,'figTitles');
addParameter(p,'relTo',2,@(x) validateattributes(x,{'numeric'},{'scalar','positive','integer'}));
addParameter(p,'normLogLike',false,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'pltYvalTxt',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'iSub',[],@(x) validateattributes(x,{'numeric'},{'scalar','positive','integer'}));

parse(p,fitRes,models,dataSubj,condVec,partVec,runEffect,iROI,modelsToPlot,...
      hFig,figTitles,varargin{:});

fitRes = p.Results.fitRes;
models = p.Results.models;
dataSubj = p.Results.dataSubj;
condVec= p.Results.condVec;
partVec = p.Results.partVec;
runEffect = p.Results.runEffect;
iROI = p.Results.iROI;
modelsToPlot = p.Results.modelsToPlot;
hFig = p.Results.hFig;
figTitles = p.Results.figTitles;
% Default is to plot relative to the first model which is not the null
relTo = p.Results.relTo;
normLogLike = p.Results.normLogLike;
pltYvalTxt = p.Results.pltYvalTxt;
iSub = p.Results.iSub;

Tgroup = fitRes{iROI,1};
Tcross = fitRes{iROI,2};
thetaCr = fitRes{iROI,4};
if size(fitRes,2) > 5
    % If the cross-validated G matrix is available, plot that
    Gpred = fitRes{iROI,6};
    isGroupFit = true;
else
    % If the cross-validated G matrix is not available, plot the
    % group-fitted one
    Gpred = fitRes{iROI,5};
    isGroupFit = false;
end

% Plotting fitting results
set(0,'currentfigure',hFig(1));
set(gcf,'Position',[40 300 1820 620]);
subplot(3,numel(figTitles),iROI);
if ~isempty(iSub)
    plotLikeRel(Tgroup,modelsToPlot,1,iSub);
else
    if isGroupFit
        T = pcm_plotModelLikelihood(Tcross,models,'Nnull',1,...
            'mindx',modelsToPlot,...
            'upperceil',Tgroup.likelihood(:,end),...
            'normalize',normLogLike);
        
        % Plot mean values above bars
        if pltYvalTxt
            Y = mean(T.likelihood_norm(:,modelsToPlot))';
            ycoord = Y;
            ycoord(ycoord < 0) = 0;
            text(1:length(Y),ycoord,num2str(Y,'%0.1f'),'HorizontalAlignment','center',...
                'VerticalAlignment','bottom','FontSize',7,'Margin',0.01,...
                'BackgroundColor','w');
        end
    else
        plotLikeRel(Tgroup,modelsToPlot,1);
    end
end
set(gca,'XTick',[]);
if iROI == numel(figTitles)
    h = findobj(gca,'Type','Bar');
    legend(flipud(h),arrayfun(@(x) models{x}.name,modelsToPlot,'UniformOutput',false));
end
title(figTitles{iROI});
if iROI > 1
    ylabel('');
end
% Plotting relative logLike relative to the required model
subplot(3,numel(figTitles),numel(figTitles)+iROI);
if ~isempty(iSub)
    plotLikeRel(Tgroup,modelsToPlot,relTo,iSub);
else
    if isGroupFit
        T = pcm_plotModelLikelihood(Tcross,models,'mindx',modelsToPlot,'Nnull',relTo,...
            'upperceil',Tgroup.likelihood(:,end),...
            'normalize',normLogLike);
        % Plot mean values above bars
        if pltYvalTxt
            Y = mean(T.likelihood_norm(:,modelsToPlot))';
            ycoord = Y;
            ycoord(ycoord < 0) = 0;
            text(1:length(Y),ycoord,num2str(Y,'%0.1f'),'HorizontalAlignment','center',...
                'VerticalAlignment','bottom','FontSize',7,'Margin',0.01,...
                'BackgroundColor','w');
        end
    else
        plotLikeRel(Tgroup,modelsToPlot,relTo);
    end
end
set(gca,'XTick',[]);
if iROI == 1
    ylabel(sprintf('logLike rel %s',models{relTo}.name));
else
    ylabel('');
end

% Plotting thetas
subplot(3,numel(figTitles),2*numel(figTitles)+iROI);
if ~isempty(iSub)
    plotThetas(thetaCr,modelsToPlot,iSub);
else
    plotThetas(thetaCr,modelsToPlot);
end
if iROI == 1
    ylabel('Thetas');
end

plotMode = {'G','RDM'};
plotFreeNonCv = false;
for i = 1:numel(plotMode)
    nRowsPlt = 2+numel(modelsToPlot)+double(plotFreeNonCv);
    % Plotting fitted models
    set(0, 'currentfigure', hFig(i+1));
    set(gcf,'Position',[300 40 1260 960]);
    % Plotting G matrix/RDM predicted from the data
    subplot(nRowsPlt,numel(figTitles),iROI);
    if ~isempty(iSub)
        plotDataG(dataSubj,partVec,condVec,runEffect,plotMode{i},iSub);
    else
        plotDataG(dataSubj,partVec,condVec,runEffect,plotMode{i});
    end
    title(figTitles{iROI});
    if iROI == 1
        ylabel(sprintf('%s\nObserved',plotMode{i}));
    end
    % Plotting predicted model G matrices/RDMs
    for j = 1:numel(modelsToPlot)
        subplot(nRowsPlt,numel(figTitles),j*numel(figTitles)+iROI);
        if ~isempty(iSub)
            plotModel(Gpred{modelsToPlot(j)},plotMode{i},iSub)
        else
            plotModel(Gpred{modelsToPlot(j)},plotMode{i})
        end
        if iROI == 1
            ylabel(sprintf('%s\n%s',plotMode{i},models{modelsToPlot(j)}.name));
        end
    end
    % Plotting freechol matrix
    subplot(nRowsPlt,numel(figTitles),(j+1)*numel(figTitles)+iROI);
    plotModel(Gpred{numel(models)},plotMode{i})
    if iROI == 1
        ylabel(sprintf('%s\n%s',plotMode{i},models{numel(models)}.name));
    end
    
    if plotFreeNonCv
        % Plotting freechol matrix non-crossvalidated
        subplot(nRowsPlt,numel(figTitles),(j+2)*numel(figTitles)+iROI);
        plotModel(fitRes{iROI,5}{numel(models)},plotMode{i})
        if iROI == 1
            ylabel(sprintf('%s\n%s\nnonCV',plotMode{i},models{numel(models)}.name));
        end
    end
end

end

function plotLikeRel(T,toPlot,relTo,varargin)
% This is a collection of snippets from pcm_plotModelLikelihood to plot
% relative likelihoods simply without all the options of that function. 
if ~isempty(varargin), iSub = varargin{1}; else, iSub = []; end
colors    = {[.7 0 0],...      % red
             [0 0 .7],...      % blue
             [.9 .6 0],...     % orange
             [0 0.6 0.6],...   % cyan
             [0.5 0 0.5],...   % purple
             [0.2 0.6 0.2],... % green
             [.7 0 0],...      % red
             [0 0 .7],...      % blue
             [.9 .6 0],...     % orange
             [0 0.6 0.6],...   % cyan
             [0.5 0 0.5],...   % purple
             [0.2 0.6 0.2]};   % green   
% ceilcolor = [0.8 0.8 0.8];
plot_fcn = @(x,y)  bar(x,y,'FaceColor',colors{x},'EdgeColor',colors{x});
ebar_fcn = @(x,ub,lb) line([x,x-0.05,x-0.05; x,x+0.05,x+0.05],[ub,ub,lb; lb,ub,lb],...
                                   'LineWidth',2,...
                                   'Color',[0 0 0]);
% Correct the upper noise ceiling for the new zero point.
% upperceil = upperceil - T.likelihood(:,relTo);
temp = T.likelihood-repmat(T.likelihood(:,relTo),1,size(T.likelihood,2));
if isempty(iSub)
    % Plot mean likelihood
    meanLike = mean(temp);
    semLike = std(temp)/sqrt(size(temp,1));
else
    % plot individual likelihood
    meanLike = temp(iSub,:);
    semLike = [];
end
% lowerceil = temp(:,end);
for j = 1:numel(toPlot)
    plot_fcn(j,meanLike(toPlot(j)));
    hold on;
    if ~isempty(semLike)
        ebar_fcn(j,meanLike(toPlot(j))+semLike(toPlot(j)),...
            meanLike(toPlot(j))-semLike(toPlot(j)));
    end
    text(j,max([meanLike(toPlot(j)),0]),sprintf('%.1f',meanLike(toPlot(j))),...
         'vert','bottom','horiz','center','FontSize',7,'Margin',0.01,...
            'BackgroundColor','w');
end
% Draw a patch encompassing area between upper and lower noise ceilings
% v = [0,mean(lowerceil); 0,mean(upperceil);...
%      numel(toPlot)+1,mean(upperceil); numel(toPlot)+1,mean(lowerceil)];
% f = [1:4];
% patch('Vertices',v,'Faces',f,'EdgeColor',ceilcolor,...
%     'FaceColor',ceilcolor,'FaceAlpha',.75);
% draw lower noise ceiling as black dotted line
% line([0;numel(toPlot)+1],[mean(lowerceil);mean(lowerceil)],...
%     'LineWidth',1.5,'Color','k','LineStyle','-.');
% Small differences are anticipated, so I set ylim relatively small
yl = ylim;
% yl = max(abs(ylim));
ylim([yl(1),1.2*yl(2)]);
set(gca,'XTick',[]);
                               
end

function plotThetas(theta,modelsToPlot,varargin)
if ~isempty(varargin), iSub = varargin{1}; else, iSub = []; end
colors    = {[.7 0 0],...      % red
             [0 0 .7],...      % blue
             [.9 .6 0],...     % orange
             [0 0.6 0.6],...   % cyan
             [0.5 0 0.5],...   % purple
             [0.2 0.6 0.2]};   % green
colors = repmat(colors,1,3);

for i = 1:numel(modelsToPlot)
    if isempty(iSub)
        m = mean(theta{modelsToPlot(i)},2);
        e = std(theta{modelsToPlot(i)},[],2)/sqrt(size(theta{modelsToPlot(i)},2));
    else
        m = theta{modelsToPlot(i)}(:,iSub);
        e = [];
    end
    if numel(m) > 1
        x = linspace(-0.4,0.4,numel(m));
    else
        x = 0;
    end
    if isempty(e)
        plot(x+i,m,'Marker','o','Color',colors{i},...
            'MarkerEdgeColor',colors{i}); hold on;
    else
        errorbar(x+i,m,e,'Marker','o','Color',colors{i},...
            'MarkerEdgeColor',colors{i}); hold on;
    end
    set(gca,'XTick',[]);
end

end

function plotDataG(Y,partVec,condVec,runEffect,plotMode,varargin)
% ----------------------------------------------------------------
% Data visualisation 
% 1. visualize the representational stucture: Get the crossvalidated
%    G-matrix for each data set and plot the average of this
if ~isempty(varargin), iSub = varargin{1}; else, iSub = []; end
if isempty(iSub)
    for s=1:length(Y)
        G_hat(:,:,s)=pcm_estGCrossval(Y{s},partVec{s},condVec{s});
    end
    Gm = mean(G_hat,3); % Mean estimate  
else
    Gm = pcm_estGCrossval(Y{iSub},partVec{iSub},condVec{iSub});
end
nConds = numel(unique(condVec{end}));
if strcmp(runEffect,'fixed')
    % Center G matrix if the runEffect is set to 'fixed'
    H = eye(nConds)-ones(nConds)/nConds;
    G = H*Gm*H';
else
    G = Gm;
end

if strcmp(plotMode,'G')
    imagesc(G);
    set(gca,'XTick',[],'YTick',[]);
    axis square;
    colorbar;
elseif strcmp(plotMode,'RDM')
    C = pcm_indicatorMatrix('allpairs',[1:nConds]');
    RDM = squareform(diag(C'*G*C));
    imagesc(RDM);
    set(gca,'XTick',[],'YTick',[]);
    axis square;
    colorbar;
end

end

function plotModel(Gpred,plotMode,varargin)
if ~isempty(varargin), iSub = varargin{1}; else, iSub = []; end
[numCond,~,numSubj] = size(Gpred);
if strcmp(plotMode,'G')
    % Averaging Gs over subjects
    if isempty(iSub)
        G = mean(Gpred,3);
    else
        G = Gpred(:,:,iSub);
    end
    imagesc(G);
    set(gca,'XTick',[],'YTick',[]);
    axis square
    colorbar;
elseif strcmp(plotMode,'RDM')
    C = pcm_indicatorMatrix('allpairs',[1:numCond]');
    if isempty(iSub)
        % Converting G to RDM separately for each subject
        temp = arrayfun(@(x) squareform(diag(C'*Gpred(:,:,x)*C)),1:numSubj,...
            'UniformOutput',false);
        % Averaging RDMs over subjects
        RDM = mean(cat(3,temp{:}),3);
    else
        RDM = squareform(diag(C'*Gpred(:,:,iSub)*C));
    end
    imagesc(RDM);
    set(gca,'XTick',[],'YTick',[]);
    axis square
    colorbar;
end
end