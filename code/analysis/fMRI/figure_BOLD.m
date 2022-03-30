function bold_table = figure_BOLD(dataPath,varargin)
% Plot mean bold responses across auditory locations and ROIs

% Parsing input
p = inputParser;
addParameter(p,'plotSingleSubjectBold',false,...
                @(x) validateattributes(x,{'logical'},{'scalar'}));
parse(p,varargin{:});
plotSingleSubjectBold = p.Results.plotSingleSubjectBold;

% Loading data
data = load(dataPath);
gimg = data.out.gimg;
aloc = data.out.aloc;
conds = data.out.conds;
roiNames = data.out.roiNames;

% Plot single subject BOLD
nconds = numel(conds);
nlocs = numel(aloc);
nsubjects = size(gimg,3);
nROIs = numel(roiNames);
col = {[133,192,249]./255,[169,90,161]./255,[15,32,128]./255};
if plotSingleSubjectBold
    for i=1:nsubjects
        figure;
        set(gcf, 'Position', [0 0 1000 200], 'PaperPositionMode', 'auto');
        for j=1:nROIs
            subplot(nROIs/5,5,j);
            for k=1:nconds
                plot(aloc,gimg(:,k,i,j),'Color',col{k}); hold on;
                set(gca, 'XTick', aloc, 'XLim', [-13 13], 'YLim', [-1.5 4]);
            end
            title(roiNames{j});
        end
    end
end
gg = squeeze(mean(gimg, 3));
gg_std = squeeze(std(gimg,[],3));
gg_sem = gg_std./(sqrt(size(gimg,3)));

% Plot group mean BOLD with sem
figure;
set(gcf, 'Position', [0 0 1000 200], 'PaperPositionMode', 'auto');
for j=1:nROIs
    subplot(nROIs/5,5,j);
    for k=1:nconds
        h{k} = shadedErrorBar(aloc,gg(:,k,j),gg_sem(:,k,j),{'Color',col{k},'LineWidth',1.5},1); hold on;
        set(gca, 'XTick', aloc, 'XLim', [-13 13],  'YLim', [-1 2.5]);
    end
    title(roiNames{j});
    if j == nROIs
        h = cat(1,h{:});
        legend([h.mainLine],conds);
    end
end

% Assign figure data to table
temp_indiv = permute(gimg,[1,2,4,3]);
s = size(temp_indiv);
temp_indiv = reshape(temp_indiv,s(1)*s(2)*s(3),s(4));
bold_table = table(reshape(repmat(roiNames,nlocs*nconds,1),[],1),...
    reshape(repmat(conds,nlocs,nROIs),[],1),...
    repmat(aloc',nconds*nROIs,1),gg(:),gg_sem(:),temp_indiv,...
    'VariableNames',{'roi','condition','auditory_location','mean_bold',...
                     'sem_bold','indiv'});
end