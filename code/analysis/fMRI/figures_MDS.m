% Loading model RDMs
dataPath = fullfile(get_path('project'),'results','data','RDM_mdl.mat');
if ~exist(dataPath,'file')
    run_generate_model_RDMs();
end
load(dataPath);

ranktransform = 'indep';

% Hemifield model MDS 
RDM_mdl_hemi_post = rsa.util.rankTransform_indepParts(RDM_mdl_hemi_post_shift,4);
MDS_mdl_hemi_post = mdscale(RDM_mdl_hemi_post,1,'Criterion','stress');

RDM_mdl_hemi_pre = rsa.util.scale01(rsa.util.rankTransform_equalsStayEqual(RDM_mdl_hemi_pre,1));
MDS_mdl_hemi_pre = mdscale(RDM_mdl_hemi_pre,1,'Criterion','stress');

% Hemifield model MDS 
RDM_mdl_dec_post = rsa.util.rankTransform_indepParts(RDM_mdl_dec_post_shift,4);
MDS_mdl_dec_post = mdscale(RDM_mdl_dec_post,1,'Criterion','stress');

RDM_mdl_dec_pre = rsa.util.scale01(rsa.util.rankTransform_equalsStayEqual(RDM_mdl_dec_pre,1));
MDS_mdl_dec_pre = mdscale(RDM_mdl_dec_pre,1,'Criterion','stress');

% Load fMRI MDS
dataPath = fullfile(get_path('project'),'results','data','RDM_fmri.mat');
if ~exist(dataPath,'file')
    RDM_fmri = get_RDMs_fMRI('iseucnorm',false);
    save(dataPath,'RDM_fmri');
end
load(dataPath,'RDM_fmri');

RDM_fMRI_post = cellfun(@(x) mean(x,3),RDM_fmri.post,'UniformOutput',false);
switch ranktransform
    case 'indep'
        % Rank transform independently
        RDM_fMRI_post = cellfun(@(x) rsa.util.rankTransform_indepParts(x,4),...
            RDM_fMRI_post,'UniformOutput',false);
    case 'all'
        % Rank transform all together
        RDM_fMRI_post = cellfun(@(x) squareform(rsa.util.rankTransform(squareform(x))),...
            RDM_fMRI_post,'UniformOutput',false);
    otherwise
        % Dont do ranktransform
end
RDM_fMRI_pre = cellfun(@(x) mean(x,3),RDM_fmri.pre,'UniformOutput',false);
RDM_fMRI_pre = cellfun(@(x) rsa.util.scale01(rsa.util.rankTransform_equalsStayEqual(x,1)),...
                   RDM_fMRI_pre,'UniformOutput',false);

rfx = false;
if rfx
    temp = cellfun(@(x) squeeze(mat2cell(x,size(x,1),size(x,2),ones(size(x,3),1)))',...
                   RDM_fmri.post,'UniformOutput',false);
    temp = cat(1,temp{:});
    MDS_fMRI_post = cellfun(@(x) mdscale(x,1),temp,'UniformOutput',false);
    MDS_fMRI_post = arrayfun(@(x) mean(cat(3,MDS_fMRI_post{x,:}),3),(1:5)', ...
                         'UniformOutput',false);
else
    MDS_fMRI_post = cellfun(@(x) mdscale(x,1,'Criterion','stress'),...
                        RDM_fMRI_post,'UniformOutput',false);
    MDS_fMRI_pre = cellfun(@(x) mdscale(x,1,'Criterion','stress'),...
                        RDM_fMRI_pre,'UniformOutput',false);
end

% Masks for extracting relevant parts of RDMs
roiNames = {'HG', 'hA', 'IPL', 'IPS', 'FEF'};
aloc = [-12,-5,-2,0,2,5,12]';
alocText = cellfun(@num2str,num2cell(repmat(aloc,2,1)),'UniformOutput',false);

% MDS plot post-test
cmap = viridis;
% cmap = magma;
colors = cmap(floor(linspace(1,256,7)),:);
figure();
nCols = size(MDS_fMRI_post,1);
xrev_post = false(1,nCols);
for i = 1:nCols
    subplot(2,nCols,i);
    imagesc(RDM_fMRI_post{i});
    set(gca,'XTick',[4,11],'XTickLabels',{'VA','AV'},...
            'YTick',[4,11],'YTickLabels',{'VA','AV'});
    axis square
    colormap(magma);
    % Computing grid line coordinates
    width = size(RDM_fMRI_post{i},1);
    x = (0:round(width/2):width)+0.5;
    x([1,end]) = [];
    x = repmat(x,2,1);
    y = repmat([0,width+0.5]',1,2);
    % Drawing vertical grid lines
    line(x,y,'Color','w','LineWidth',1.5);
    % Drawing horizontal grid lines
    line(y,x,'Color','w','LineWidth',1.5);
    title(roiNames{i});
    if i == nCols
        c = colorbar();
        c.Label.String = 'Normalised distance';
    end
    
    subplot(2,nCols,i+nCols);
    scatter(MDS_fMRI_post{i}(1:7),0.2*ones(7,1),96,colors,'filled',...
            'MarkerEdgeColor','k'); hold on;
    scatter(MDS_fMRI_post{i}(8:14),-0.2*zeros(7,1),96,colors,'filled',...
            'MarkerEdgeColor','k');
    text([mean(MDS_fMRI_post{i}(1:7)),mean(MDS_fMRI_post{i}(8:14))],[0.35;-0.28],...
        {'VA-adaptation','AV-adaptation'},'HorizontalAlignment','center',...
        'VerticalAlignment','bottom');
    if i == nCols
        colormap(gca,'viridis');
        c = colorbar('Ticks',floor(linspace(1,256,7))/256, ...
                     'TickLabels',alocText,'Location','manual','Position',...
                     [0.95,0.1,0.02,0.3]);
        c.Label.String = 'Auditory location (deg)';
    end
    if MDS_fMRI_post{i}(1) > MDS_fMRI_post{i}(7)
        xrev_post(i) = true;
        set(gca,'XDir','reverse');
    end
    % title(labels{i});
    set(gca,'Ytick',[]);
    ylim([-0.6,0.6]);
    if isempty(ranktransform) || strcmp(ranktransform,'all')
        if xrev_post(i)
            xlims = [min(MDS_fMRI_post{i}(8:14)),max(MDS_fMRI_post{i}(1:7))];
            breaklims = [max(MDS_fMRI_post{i}(8:14)),min(MDS_fMRI_post{i}(1:7))];
        else
            xlims = [min(MDS_fMRI_post{i}(1:7)),max(MDS_fMRI_post{i}(8:14))];
            breaklims = [max(MDS_fMRI_post{i}(1:7)),min(MDS_fMRI_post{i}(8:14))];
        end
        if xlims(1) > xlims(2), xlims = fliplr(xlims); end
        if breaklims(1) > breaklims(2), breaklims = fliplr(breaklims); end
        xlim(xlims);
        % Break x axis for better visibility
        breakxaxis(breaklims,0.0025);
    end
    
end

% Hemifield model MDS post-test
figure();
subplot(2,1,1);
imagesc(RDM_mdl_hemi_post);
set(gca,'XTick',[4,11],'XTickLabels',{'VA','AV'},...
        'YTick',[4,11],'YTickLabels',{'VA','AV'});
axis square
% Computing grid line coordinates
width = size(RDM_mdl_hemi_post,1);
x = (0:round(width/2):width)+0.5;
x([1,end]) = [];
x = repmat(x,2,1);
y = repmat([0,width+0.5]',1,2);
% Drawing vertical grid lines
line(x,y,'Color','w','LineWidth',1.5);
% Drawing horizontal grid lines
line(y,x,'Color','w','LineWidth',1.5);
title('Hemifield model');
c = colorbar();
c.Label.String = 'Normalised distance';

subplot(2,1,2);
scatter(MDS_mdl_hemi_post(1:7),0.2*ones(7,1),96,colors,'filled',...
        'MarkerEdgeColor','k'); hold on;
scatter(MDS_mdl_hemi_post(8:14),-0.2*zeros(7,1),96,colors,'filled',...
        'MarkerEdgeColor','k');
text([0,0],[0.35;-0.28],{'VA-adaptation','AV-adaptation'},...
     'HorizontalAlignment','center','VerticalAlignment','bottom');
colormap(gca,'viridis');
c = colorbar('Ticks',floor(linspace(1,256,7))/256, ...
             'TickLabels',alocText);
c.Label.String = 'Auditory location (deg)';

if MDS_mdl_hemi_post(1) > MDS_mdl_hemi_post(7)
    set(gca,'XDir','reverse')
end   
% title(labels{i});
set(gca,'Ytick',[]);
ylim([-0.6,0.6]);

% Decisional model MDS post-test
figure();
subplot(2,1,1);
imagesc(RDM_mdl_dec_post);
set(gca,'XTick',[4,11],'XTickLabels',{'VA','AV'},...
        'YTick',[4,11],'YTickLabels',{'VA','AV'});
axis square
% Computing grid line coordinates
width = size(RDM_mdl_dec_post,1);
x = (0:round(width/2):width)+0.5;
x([1,end]) = [];
x = repmat(x,2,1);
y = repmat([0,width+0.5]',1,2);
% Drawing vertical grid lines
line(x,y,'Color','w','LineWidth',1.5);
% Drawing horizontal grid lines
line(y,x,'Color','w','LineWidth',1.5);
title('Decisional model');
c = colorbar();
c.Label.String = 'Normalised distance';

subplot(2,1,2);
scatter(MDS_mdl_dec_post(1:7),0.2*ones(7,1),96,colors,'filled',...
        'MarkerEdgeColor','k'); hold on;
scatter(MDS_mdl_dec_post(8:14),-0.2*zeros(7,1),96,colors,'filled',...
        'MarkerEdgeColor','k');
text([0,0],[0.35;-0.28],{'VA-adaptation','AV-adaptation'},...
     'HorizontalAlignment','center','VerticalAlignment','bottom');
colormap(gca,'viridis');
c = colorbar('Ticks',floor(linspace(1,256,7))/256, ...
             'TickLabels',alocText);
c.Label.String = 'Auditory location (deg)';

% if MDS_mdl_dec_post(1) > MDS_mdl_dec_post(7)
%     set(gca,'XDir','reverse')
% end   
% title(labels{i});
set(gca,'Ytick',[]);
ylim([-0.6,0.6]);

% MDS plot pre-test
cmap = viridis;
% cmap = magma;
colors = cmap(floor(linspace(1,256,7)),:);
figure();
nCols = size(MDS_fMRI_pre,1);
xrev_pre = false(1,nCols);
for i = 1:nCols
    subplot(2,nCols,i);
    imagesc(RDM_fMRI_pre{i});
    set(gca,'XTick',1:7,'XTickLabels',alocText,...
            'YTick',1:7,'YTickLabels',alocText);
    axis square
    colormap(magma);
    title(roiNames{i});
    if i == nCols
        c = colorbar();
        c.Label.String = 'Normalised distance';
    end
    
    subplot(2,nCols,i+nCols);
    scatter(MDS_fMRI_pre{i},zeros(7,1),96,colors,'filled',...
            'MarkerEdgeColor','k'); hold on;
    if i == nCols
        colormap(gca,'viridis');
        c = colorbar('Ticks',floor(linspace(1,256,7))/256, ...
                     'TickLabels',alocText);
        c.Label.String = 'Auditory location (deg)';
    end
    if MDS_fMRI_pre{i}(1) > MDS_fMRI_pre{i}(7)
        xrev_pre(i) = true;
        set(gca,'XDir','reverse')
    end   
    % title(labels{i});
    set(gca,'Ytick',[]);
    ylim([-0.6,0.6]);
end

% Computing the rank-correlation between spatial locations and the MDS
% results to see how congruent is the ordering of the MDS with the real
% locations. 
% Pre-test
% I take the absolute value of the correlation coefficient as the MDS
% sometimes flips the order of the locations. 
rho_pre = cellfun(@(c) corr(aloc,c,'Type','Spearman'),MDS_fMRI_pre);
% Post-test
rho_post = cellfun(@(c) corr(repmat(aloc,2,1),c,'Type','Spearman'),MDS_fMRI_post);
rho_MDS_aloc = table(cat(1,repmat({'pre-test'},5,1),repmat({'post-test'},5,1)),...
                     repmat(roiNames',2,1),...
                     cat(1,abs(rho_pre),abs(rho_post)),...
                     'VariableNames',{'Exp_phase','ROI','Rho'});

% Hemifield model pre-test
figure();
subplot(2,1,1);
imagesc(RDM_mdl_hemi_pre);
set(gca,'XTick',1:7,'XTickLabels',alocText,...
            'YTick',1:7,'YTickLabels',alocText)
axis square
title('Hemifield model');
c = colorbar();
c.Label.String = 'Normalised distance';

subplot(2,1,2);
scatter(MDS_mdl_hemi_pre,zeros(7,1),96,colors,'filled',...
        'MarkerEdgeColor','k'); hold on;
colormap(gca,'viridis');
c = colorbar('Ticks',floor(linspace(1,256,7))/256, ...
             'TickLabels',alocText);
c.Label.String = 'Auditory location (deg)';

if MDS_mdl_hemi_pre(1) > MDS_mdl_hemi_pre(7)
    set(gca,'XDir','reverse')
end   
% title(labels{i});
set(gca,'Ytick',[]);
ylim([-0.6,0.6]);

% Decisional model pre-test
figure();
subplot(2,1,1);
imagesc(RDM_mdl_dec_pre);
set(gca,'XTick',1:7,'XTickLabels',alocText,...
            'YTick',1:7,'YTickLabels',alocText)
axis square
title('Decisional model');
c = colorbar();
c.Label.String = 'Normalised distance';

subplot(2,1,2);
scatter(MDS_mdl_dec_pre,zeros(7,1),96,colors,'filled',...
        'MarkerEdgeColor','k'); hold on;
colormap(gca,'viridis');
c = colorbar('Ticks',floor(linspace(1,256,7))/256, ...
             'TickLabels',alocText);
c.Label.String = 'Auditory location (deg)';

% if MDS_mdl_dec_pre(1) > MDS_mdl_dec_pre(7)
%     set(gca,'XDir','reverse')
% end   
% title(labels{i});
set(gca,'Ytick',[]);
ylim([-0.6,0.6]);

%% Saving MDS plot data to table
temp_roi = repmat(roiNames,numel(aloc),1);
temp_roi = temp_roi(:);
temp_cond = repmat({'pre'},numel(temp_roi),1);
temp_aloc = repmat(aloc,numel(roiNames),1);
temp_mds = cat(1,MDS_fMRI_pre{:});
mds_pre_table = table(temp_roi,temp_cond,temp_aloc,temp_mds,'VariableNames',...
                      {'roi','cond','aloc','mds'});
temp_roi = repmat(roiNames,numel(aloc)*2,1);
temp_roi = temp_roi(:);
temp_cond = repmat({'postVA','postAV'},numel(aloc),1);
temp_cond = repmat(temp_cond(:),numel(roiNames),1);
temp_aloc = repmat(cat(1,aloc,aloc),numel(roiNames),1);
temp_mds = cat(1,MDS_fMRI_post{:});
mds_post_table = table(temp_roi,temp_cond,temp_aloc,temp_mds,'VariableNames',...
                      {'roi','cond','aloc','mds'});
mds_table = cat(1,mds_pre_table,mds_post_table);