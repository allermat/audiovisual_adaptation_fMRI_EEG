function logBF_table = figure_lme_fMRI_pre(dataPath,varargin)

% Parsing input
p = inputParser;

validThirdModels = {'Pmf','Rt'};
validHemispheres = {'lh','rh','mean_lhrh','lh_and_rh'};

addRequired(p,'dataPath',@ischar);
addParameter(p,'hemisphere','lh_and_rh',@(x) ismember(x,validHemispheres));
addParameter(p,'thirdModel','',@(x) ismember(x,validThirdModels));
addParameter(p,'pltYvalTxt',true,@(x) validateattributes(x,{'logical'},{'scalar'}));

parse(p,dataPath,varargin{:});

dataPath = p.Results.dataPath;
thirdModel = p.Results.thirdModel;
pltYvalTxt = p.Results.pltYvalTxt;

data = load(dataPath);
BIC_rel = data.out.BIC_rel;
roiNames = BIC_rel.Properties.VariableNames;
nROIs = numel(roiNames);
modelNames = BIC_rel.Properties.RowNames;

figure();
set(gcf, 'Position', [0 0 1000 200], 'PaperPositionMode', 'auto');

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
% Set y-axis limits
temp = BIC_rel{:,:};
yl = [min(temp(:)),max(temp(:))]*1.05;

% Which models to include in the plot
if isempty(thirdModel)
    plotOrder = 2:4;
elseif strcmp(thirdModel,'Pmf')
    plotOrder = 2:8;
else
    error('Not yet implemented!');
end

for iROI = 1:nROIs
    subplot(1,nROIs,iROI);
    for k = 1:numel(plotOrder)
        bar(k,BIC_rel{plotOrder(k),iROI},'FaceColor',colors{k}); hold on;
    end
    if pltYvalTxt
        Y = BIC_rel{plotOrder,iROI};
        ycoord = Y;
        ycoord(ycoord < 0) = 0;
        text(1:length(Y),ycoord,num2str(Y,'%0.1f'),'HorizontalAlignment','center',...
            'VerticalAlignment','bottom','FontSize',7,'Margin',0.01,...
            'BackgroundColor','w');
    end
    if iROI == 1, ylabel('logBF'); end
    ylim(yl);
    set(gca,'XTickLabels',[]);
    title(roiNames{iROI});
end
legend(BIC_rel.Properties.RowNames(plotOrder));

temp_roi = repmat(roiNames,numel(plotOrder),1);
temp_roi = temp_roi(:);
temp_model = repmat(modelNames(plotOrder),numel(roiNames),1);
temp_data = BIC_rel{plotOrder,:};
logBF_table = table(temp_roi,temp_model,temp_data(:),'VariableNames',...
                    {'roi','model','logBF'});
end