function varargout = plot_group_RT(varargin)
% Function for plotting response time data 
% 
% USAGE: 
%   plot_group_RT()
%   plot_group_RT(Name,Value) 
%
% INPUT (Name-Value pair arguments): 
%   plotGroupError: whether to plot group-level errors and in which form.
%       Possible values are 'shaded','errorbar','none'. Default: 'shaded'.
%   removeMeanAcrossLoc: whether to subtract mean across locations for
%       group-level RTs. Default: true. 
%
% OUTPUT: 
%   Response time tables, one for each modality: psychophysics, fMRI, EEG
%   For each modality, it draws one group level and one individual level
%       RT figure. 
% 
% Mate Aller, 2021, allermat@gmail.com

validErrors = {'shaded','errorbar','none'};

% Parsing input
p = inputParser;

addParameter(p,'plotGroupError','shaded',@(x) ismember(x,validErrors));
addParameter(p,'removeMeanAcrossLoc',true,...
             @(x) validateattributes(x,{'logical'},{'scalar'}));

parse(p,varargin{:});

plotGroupError = p.Results.plotGroupError;
removeMeanAcrossLoc = p.Results.removeMeanAcrossLoc;

modality = {'PP' 'fMRI' 'EEG'};
color = {[133,192,249]./255,[169,90,161]./255,[15,32,128]./255};
[indiv_RT_pre_cell,indiv_RT_post_cell,...
    group_RT_pre_cell,group_RT_post_cell] = deal(cell(1,numel(modality)));
for m = 1:numel(modality)
    % Load data
    load([modality{m} '_behav_all.mat']);
    nSub = numel(out);
    pre = cat(2, out.pre);
    post = cat(2, out.post);
    data = table();
    for iSub = 1:nSub
        % In 1 subject there are empty cells in the data because there was
        % only 2 sesson of each adaptation direction, getting rid of empty
        % cells. 
        pre(iSub).data(cellfun(@isempty,pre(iSub).data)) = [];
        post(iSub).data(cellfun(@isempty,post(iSub).data)) = [];
        temp_pre_r = dataset2table(cat(1,pre(iSub).data{:,1}));
        temp_pre_r.adaptDir = repmat({'R'},size(temp_pre_r,1),1);
        temp_pre_l = dataset2table(cat(1,pre(iSub).data{:,2}));
        temp_pre_l.adaptDir = repmat({'L'},size(temp_pre_l,1),1);
        temp_post_r = dataset2table(cat(1,post(iSub).data{:,1}));
        temp_post_r.adaptDir = repmat({'R'},size(temp_post_r,1),1);
        temp_post_l = dataset2table(cat(1,post(iSub).data{:,2}));
        temp_post_l.adaptDir = repmat({'L'},size(temp_post_l,1),1);
        temp = cat(1,temp_pre_r,temp_pre_l,temp_post_r,temp_post_l);
        temp.subID = iSub*ones(size(temp,1),1);
        data = cat(1,data,temp);
    end
    RT_pre = varfun(@nanmean,...
                    data(ismember(data.blocktype,'pre-test') & data.catch_trial == 1,:),...
                    'InputVariables',{'RT'},'GroupingVariables',...
                    {'subID','aloc'});
    RT_pre.Properties.VariableNames{'nanmean_RT'} = 'RT';
    RT_post = varfun(@nanmean,...
                     data(ismember(data.blocktype,'post-test') & data.catch_trial == 1,:),...
                     'InputVariables',{'RT'},'GroupingVariables',...
                     {'subID','adaptDir','aloc'});
    RT_post.Properties.VariableNames{'nanmean_RT'} = 'RT';
    % Plot subject specific RTs
    figure('Name',sprintf('%s_subject_level',modality{m}));
    if nSub > 5
        set(gcf, 'PaperPositionMode', 'auto', 'Position', [300 0 1600 900]);
    else
        set(gcf, 'PaperPositionMode', 'auto', 'Position', [300 300 1600 300]);
    end
    for iSub = 1:nSub
        if nSub > 5
            subplot(3,nSub/3,iSub);
        else
            subplot(1,nSub,iSub);
        end
        actSubIdx_pre = RT_pre.subID == iSub;
        actSubIdx_post = RT_post.subID == iSub;
        isRight = ismember(RT_post.adaptDir,'R');
        isLeft = ismember(RT_post.adaptDir,'L');
        plot(RT_pre.aloc(actSubIdx_pre),RT_pre.RT(actSubIdx_pre),...
             'Color',color{1},'LineWidth',1.5); hold on
        plot(RT_post.aloc(actSubIdx_post & isRight),...
             RT_post.RT(actSubIdx_post & isRight),'Color',color{2},...
             'LineWidth',1.5);
        plot(RT_post.aloc(actSubIdx_post & isLeft),...
             RT_post.RT(actSubIdx_post & isLeft),'Color',color{3},...
             'LineWidth',1.5);
        StimLevels = unique(RT_pre.aloc);
        set(gca, 'FontSize', 12, 'Xtick', StimLevels, 'box', 'off');
        xlabel(gca, 'A location (°)', 'FontSize', 16);
        ylabel(gca, 'Response time (s)', 'FontSize', 16);
        if exist('figname', 'var')
            title(figname, 'FontSize', 18);
        end
        title(sprintf('sub-%02d',iSub));
        xlim(gca, [-15 15]);
        ylim(gca, [min(cat(1,RT_pre.RT,RT_post.RT)),...
                   max(cat(1,RT_pre.RT,RT_post.RT))]);
        if iSub == nSub
            legend({'pre-adaptation','postAV-adaptation','postVA-adaptation'},...
                   'Location','NorthEast');
        end
    end
    
    if removeMeanAcrossLoc
        % Remove mean across locations for each subject before computing group
        % average
        subj_mean_pre = varfun(@mean,RT_pre,'InputVariables','RT',...
            'GroupingVariables',{'subID'});
        subj_mean_post = varfun(@mean,RT_post,'InputVariables','RT',...
            'GroupingVariables',{'subID','adaptDir'});
        for iSub = 1:nSub
            temp = RT_pre.RT(RT_pre.subID == iSub);
            RT_pre.RT(RT_pre.subID == iSub) = ...
                temp-subj_mean_pre.mean_RT(subj_mean_pre.subID == iSub);
            
            actSel_RT = RT_post.subID == iSub & strcmp(RT_post.adaptDir,'L');
            actSel_meanRT = subj_mean_post.subID == iSub & ...
                            strcmp(subj_mean_post.adaptDir,'L');
            RT_post.RT(actSel_RT) = ...
                RT_post.RT(actSel_RT)-subj_mean_post.mean_RT(actSel_meanRT);
            
            actSel_RT = RT_post.subID == iSub & strcmp(RT_post.adaptDir,'R');
            actSel_meanRT = subj_mean_post.subID == iSub & ...
                            strcmp(subj_mean_post.adaptDir,'R');
            RT_post.RT(actSel_RT) = ...
                RT_post.RT(actSel_RT)-subj_mean_post.mean_RT(actSel_meanRT);
            
        end
    end
    figure('Name',sprintf('%s_group_mean',modality{m}));
    group_RT_pre = join(varfun(@mean,RT_pre,'InputVariables','RT',...
                          'GroupingVariables',{'aloc'}),...
                        varfun(@std,RT_pre,'InputVariables','RT',...
                          'GroupingVariables',{'aloc'}));
    group_RT_pre.sem_RT = group_RT_pre.std_RT./sqrt(group_RT_pre.GroupCount);
    group_RT_post = join(varfun(@mean,RT_post,'InputVariables','RT',...
                           'GroupingVariables',{'adaptDir','aloc'}),...
                         varfun(@std,RT_post,'InputVariables','RT',...
                           'GroupingVariables',{'adaptDir','aloc'}));
    group_RT_post.sem_RT = group_RT_post.std_RT./sqrt(group_RT_post.GroupCount);
    isRight = ismember(group_RT_post.adaptDir,'R');
    isLeft = ismember(group_RT_post.adaptDir,'L');
    switch plotGroupError
        case 'shaded'
            h1 = shadedErrorBar(group_RT_pre.aloc,group_RT_pre.mean_RT,...
                    group_RT_pre.sem_RT,{'Color',color{1},'LineWidth',2},1); hold on
            h2 = shadedErrorBar(group_RT_post.aloc(isRight),group_RT_post.mean_RT(isRight),...
                    group_RT_post.sem_RT(isRight),{'Color',color{2},'LineWidth',2},1);
            h3 = shadedErrorBar(group_RT_post.aloc(isLeft),group_RT_post.mean_RT(isLeft),...
                    group_RT_post.sem_RT(isLeft),{'Color',color{3},'LineWidth',2},1);
            legend([h1.mainLine,h2.mainLine,h3.mainLine],...
                   {'pre-adaptation','postAV-adaptation','postVA-adaptation'},...
                   'Location','SouthEast');
        case 'errorbar'
            errorbar(group_RT_pre.aloc,group_RT_pre.mean_RT,...
                group_RT_pre.sem_RT,'Color',color{1},'LineWidth',2); hold on
            errorbar(group_RT_post.aloc(isRight),group_RT_post.mean_RT(isRight),...
                group_RT_post.sem_RT(isRight),'Color',color{2},'LineWidth',2);
            errorbar(group_RT_post.aloc(isLeft),group_RT_post.mean_RT(isLeft),...
                group_RT_post.sem_RT(isLeft),'Color',color{3},'LineWidth',2);
        case 'none'
            plot(group_RT_pre.aloc,group_RT_pre.mean_RT,'Color',color{1},...
                'LineWidth',2); hold on
            plot(group_RT_post.aloc(isRight),group_RT_post.mean_RT(isRight),...
                'LineWidth',2,'Color',color{2});
            plot(group_RT_post.aloc(isLeft),group_RT_post.mean_RT(isLeft),...
                'LineWidth',2,'Color',color{3});
            legend({'pre-adaptation','postAV-adaptation','postVA-adaptation'},...
                   'Location','SouthEast');
    end
    title(sprintf('Group mean response times\n%s',modality{m}));
    StimLevels = unique(RT_pre.aloc);
    set(gca, 'FontSize', 12, 'Xtick', StimLevels, 'box', 'off');
    xlabel(gca, 'A location (°)', 'FontSize', 16);
    ylabel(gca, 'Response time (s)', 'FontSize', 16);
    xlim(gca, [-15 15]);
    ylim([-0.06,0.08]);
    indiv_RT_pre_cell{m} = RT_pre;
    indiv_RT_post_cell{m} = RT_post;
    group_RT_pre_cell{m} = group_RT_pre;
    group_RT_post_cell{m} = group_RT_post;
    % Write group RT data to table
    group_RT_pre.cond = repmat({'pre'},size(group_RT_pre.aloc));
    group_RT_pre.modality = repmat(lower(modality(m)),size(group_RT_pre.aloc));
    group_RT_pre = group_RT_pre(:,[7,6,1:5]);
    temp = varfun(@(x) x',RT_pre,'InputVariables',{'RT'},'GroupingVariables',{'aloc'});
    temp = splitvars(temp,'Fun_RT','NewVariableNames',...
        strcat('sub',cellfun(@num2str,num2cell(1:size(temp.Fun_RT,2)),'UniformOutput',false),...
               '_RT'));
    group_RT_pre = join(group_RT_pre,temp);
    group_RT_post.adaptDir(ismember(group_RT_post.adaptDir,'L')) = {'postVA'};
    group_RT_post.adaptDir(ismember(group_RT_post.adaptDir,'R')) = {'postAV'};
    group_RT_post.Properties.VariableNames{1} = 'cond';
    group_RT_post.modality = repmat(lower(modality(m)),size(group_RT_post.aloc));
    group_RT_post = group_RT_post(:,[7,1:6]);
    temp = varfun(@(x) x',RT_post,'InputVariables',{'RT'},...
                  'GroupingVariables',{'adaptDir','aloc'});
    temp = splitvars(temp,'Fun_RT','NewVariableNames',...
        strcat('sub',cellfun(@num2str,num2cell(1:size(temp.Fun_RT,2)),'UniformOutput',false),...
               '_RT'));
    temp.adaptDir(ismember(temp.adaptDir,'L')) = {'postVA'};
    temp.adaptDir(ismember(temp.adaptDir,'R')) = {'postAV'};
    temp.Properties.VariableNames{1} = 'cond';
    group_RT_post = join(group_RT_post,temp);
    group_RT = cat(1,group_RT_pre,group_RT_post);
    group_RT.GroupCount = [];
    varargout{m} = group_RT;
end
