function logBF_table = figure_lme_fMRI_custom(dataPath,varargin)

p = inputParser;
addRequired(p,'dataPath',@ischar);
addParameter(p,'modelColors',{},@(x) validateattributes(x,{'cell'},{'size',[1,5]}));
addParameter(p,'pltYvalTxt',true,@(x) validateattributes(x,{'logical'},{'scalar'}));

parse(p,dataPath,varargin{:});

dataPath = p.Results.dataPath;
modelColors = p.Results.modelColors;
pltYvalTxt = p.Results.pltYvalTxt;

data = load(dataPath);
data = data.out;
BIC_rel = data.BIC_rel;
roiNames = BIC_rel.Properties.VariableNames;
nROIs = numel(roiNames);
modelNames = data.modelNames;
% plot BIC
colors    = {[.7 0 0],...      % red
    [0 0 .7],...      % blue
    [.9 .6 0],...     % orange
    [0 0.6 0.6],...   % cyan
    [0.5 0 0.5],...   % purple
    [0.2 0.6 0.2]};   % green

% Setting plotOrder to the last two models in each ROI
temp = sum(~isnan(BIC_rel{:,:}));
plotOrder = arrayfun(@(x) {x-1:x},temp);
temp = BIC_rel{:,:};
yl = [min(temp(:)),max(temp(:))]*1.05;
nFigs = size(plotOrder,1);
for iFig = 1:nFigs
    figure();
    set(gcf, 'Position', [0 0 1000 200], 'PaperPositionMode', 'auto');
    for iROI = 1:nROIs
        subplot(nFigs,nROIs,(iFig-1)*nROIs+iROI);
        for k = 1:numel(plotOrder{iFig,iROI})
            if ~isempty(modelColors)
                actColor = modelColors{iFig,iROI}{k};
            else
                actColor = colors{1};
            end
            bar(k,BIC_rel{plotOrder{iFig,iROI}(k),iROI},'FaceColor',actColor); hold on;
        end
        if pltYvalTxt
            Y = BIC_rel{plotOrder{iFig,iROI},iROI};
            ycoord = Y;
            ycoord(ycoord < 0) = 0;
            text(1:length(Y),ycoord,num2str(Y,'%0.1f'),'HorizontalAlignment','center',...
                'VerticalAlignment','bottom','FontSize',7,'Margin',0.01,...
                'BackgroundColor','w');
        end
        if iROI == 1, ylabel('logBF'); end
        ylim(yl);
        set(gca,'XTick',1:k,'XTickLabel',modelNames{iROI}(plotOrder{iFig,iROI}),...
            'XTickLabelRotation',45);
        if iFig == 1, title(roiNames{iROI}); end
    end
end
temp_roi = cellfun(@(c1,c2) repmat({c1},numel(c2),1),roiNames,plotOrder,...
                   'UniformOutput',false);
temp_roi = cat(1,temp_roi{:});
temp_model = cellfun(@(c1,c2) c1(c2)',modelNames,plotOrder,...
                     'UniformOutput',false);
temp_model = cat(1,temp_model{:});
temp_data = BIC_rel{:,:};
temp_data = arrayfun(@(x) temp_data(plotOrder{x},x),1:size(temp_data,2),...
                     'UniformOutput',false);
temp_data = cat(1,temp_data{:});
logBF_table = table(temp_roi,temp_model,temp_data(:),'VariableNames',...
                    {'roi','model','logBF'});

end