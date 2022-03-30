function varargout = figure_model_distributions(varargin)
% Function for generating model predictions for hemifield, place and
% decisional code models. 
% 
% USAGE: 
%   figure_model_distributions()
%   figure_model_distributions(Name,Value)
% 
% DETAILS: 
%   The predictions are generated for a hypothetical brain region in the
%   left hemisphere, so ipsilateral means left and  contralateral mans 
%   right. 
%
% INPUT (Name-Value pair arguments): 
%   'hemiPred': prediction for the hemifield model, can be linear ('lin') 
%       or constant ('const'). If constant, the proportion of right and
%       left tuned neurons is set 50-50%, so the mean activity predictions
%       are constant across spatial locations. If linear the proportions
%       are set to 70-30% respectively for right and left tuned neurons. 
%   'placePred': prediction for the place model, can be linear ('lin') 
%       or constant ('const'). If constant, the maxima of neural tuning 
%       curves are distributed uniformly at inervals of 1 deg. If linear, 
%       the maxima of neural tuning curves are sampled from a normal
%       distribution with mean = 90 deg, SD = 80 deg. 
%   'aloc': vector of auditory locations where the model predictions are
%       calculated. Default [-12,-5,-2,0,2,5,12]
%   'plotCorr': whether to plot correlations across model predictions.
%       Default: false. 
% 
% OUTPUT: 
%   Draws a figure which shows for each model (1)the neural tuning 
%       functions, (2) mean predicted activations across neurons and the
%       second moment of predicted activations (as a representational
%       dissimilarity matrix) across neurons with (3) and without (4)
%       recalibration. 
% 
% Mate Aller, 2021, allermat@gmail.com

preds = {'const','lin'};

% Parsing input
p = inputParser;
addParameter(p,'hemiPred','lin',@(x) ismember(x,preds));
addParameter(p,'placePred','const',@(x) ismember(x,preds));
addParameter(p,'aloc',[-12,-5,-2,0,2,5,12],...
             @(x) validateattributes(x,{'double'},{'row','increasing'}));
addParameter(p,'plotCorr',false,@(x) validateattributes(x,{'logical'},{'scalar'}));

parse(p,varargin{:});

hemiPred = p.Results.hemiPred;
placePred = p.Results.placePred;
aloc = p.Results.aloc;
plotCorr = p.Results.plotCorr;

% Hemifield model features
conds = {'pre','postVA','postAV'};
shift = 2.3; % Derived from behavioural data
sigmaHemi = 64; % Based on Salminen et al. 2009
popCenter = 90;
nNeurons = 360;
if strcmp(hemiPred,'const')
    propR = 0.5;
else
    propR = 0.7;
end
propL = 1-propR;
spanHemi = 20;
% Mean of tuning curves uniformly distributed around +90 and -90 degrees
muHemi_L = rand(round(nNeurons*propL),1)*spanHemi-(spanHemi/2)-popCenter;
muHemi_R = rand(round(nNeurons*propR),1)*spanHemi-(spanHemi/2)+popCenter;

funHemi = @(m,s) normpdf(aloc,m+s,sigmaHemi)/normpdf(m,m,sigmaHemi);
% Pretest
y = cat(1,arrayfun(funHemi,muHemi_L,zeros(size(muHemi_L)),'UniformOutput',false),...
    arrayfun(funHemi,muHemi_R,zeros(size(muHemi_R)),'UniformOutput',false));
featHemiPre = cat(1,y{:});
% VA adaptation
y = cat(1,arrayfun(funHemi,muHemi_L,shift*ones(size(muHemi_L)),'UniformOutput',false),...
    arrayfun(funHemi,muHemi_R,shift*ones(size(muHemi_R)),'UniformOutput',false));
featHemiPostVA = cat(1,y{:});
% AV adaptation
y = cat(1,arrayfun(funHemi,muHemi_L,-shift*ones(size(muHemi_L)),'UniformOutput',false),...
    arrayfun(funHemi,muHemi_R,-shift*ones(size(muHemi_R)),'UniformOutput',false));
featHemiPostAV = cat(1,y{:});

% Place code model features
sigmaPlace = 26; % Based on Salminen et al. 2009
if strcmp(placePred,'const')
    % The means of the individual neural/voxel tuning functions are sampled
    % unifromly at equal distance
    muPlace = -179:180;
else
    % The means of the individual neural/voxel tuning functions are sampled
    % from a normal distribution with mean and sd given in the input
    meanOfMu = 90;
    sdOfMu = 80;
    for i = 1:nNeurons
        m = (randn*sdOfMu)+meanOfMu;
        while m <= -180 || m > 180
            m = (randn*sdOfMu)+meanOfMu;
        end
        muPlace(i) = m;
    end
end
funPlace = @(m) normpdf(aloc,m,sigmaPlace)/normpdf(m,m,sigmaPlace);
featPlacePre = arrayfun(funPlace,muPlace,'UniformOutput',false);
featPlacePre = cat(1,featPlacePre{:});
featPlacePostVA = arrayfun(funPlace,muPlace+shift,'UniformOutput',false);
featPlacePostVA = cat(1,featPlacePostVA{:});
featPlacePostAV = arrayfun(funPlace,muPlace-shift,'UniformOutput',false);
featPlacePostAV = cat(1,featPlacePostAV{:});

% Decision model features
spanDec = 2;
sigmaDec = 10;
muDec = rand(nNeurons,1)*spanDec-(spanDec/2);
funDec = @(m,s) abs(abs(normcdf(aloc,m+s,sigmaDec)-0.5)-0.5)*2;
y = arrayfun(funDec,muDec,zeros(size(muDec)),'UniformOutput',false);
featDecPre = cat(1,y{:});
y = arrayfun(funDec,muDec,shift*ones(size(muDec)),'UniformOutput',false);
featDecPostVA = cat(1,y{:});
y = arrayfun(funDec,muDec,-shift*ones(size(muDec)),'UniformOutput',false);
featDecPostAV = cat(1,y{:});

% Decisional choice (PMF) model features
% Load PMF data
fLoaded = load(fullfile(get_path('project'),'results','data','fMRI_behav_all.mat'));
pre = cat(2, fLoaded.out.pre);
post = cat(2, fLoaded.out.post);
% Computing group average PMF values
y = arrayfun(@(i) mean(pre(i).NumPos./pre(i).OutOfNum), 1:size(pre,2),...
    'UniformOutput', false);
featPmfPre = mean(cat(1,y{:}));
% post-AV adaptation
y = arrayfun(@(i) mean(post(i).NumPos(1:2,:)./post(i).OutOfNum(1:2,:)),...
    1:size(post,2),'UniformOutput', false);
featPmfPostAV = mean(cat(1,y{:}));
% post-VA adaptation
y = arrayfun(@(i) mean(post(i).NumPos(3:4,:)./post(i).OutOfNum(3:4,:)),...
    1:size(post,2),'UniformOutput', false);
featPmfPostVA = mean(cat(1,y{:}));

% Reaction time (RT) model features
nSub = size(pre,2);
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

group_RT_pre = varfun(@mean,RT_pre,'InputVariables','RT',...
    'GroupingVariables',{'aloc'});
group_RT_pre.Properties.VariableNames{'mean_RT'} = 'RT';
group_RT_post = varfun(@mean,RT_post,'InputVariables','RT',...
    'GroupingVariables',{'adaptDir','aloc'});
group_RT_post.Properties.VariableNames{'mean_RT'} = 'RT';
featRtPre = group_RT_pre.RT';
featRtPostAV = group_RT_post.RT(ismember(group_RT_post.adaptDir,'R'))';
featRtPostVA = group_RT_post.RT(ismember(group_RT_post.adaptDir,'L'))';

nModels = 5;
nPlotsPerModel = 4;
figure('Name','Model predictions');
ax = subplot(nPlotsPerModel,nModels,1);
set(gca, 'FontName', 'Arial');
[x_hemi,y_hemi] = plotHemi(ax,propR);
ylabel('Relative neural activity');
title('Hemifield code');

ax = subplot(nPlotsPerModel,nModels,2);
set(gca, 'FontName', 'Arial');
[x_place,y_place] = plotPlace(ax,muPlace(1:15:numel(muPlace)));
title('Place code');

ax = subplot(nPlotsPerModel,nModels,3);
set(gca, 'FontName', 'Arial');
[x_dec,y_dec] = plotDec(ax);
title('Decisional uncertainty code');

% There is no neural data for the PMF and RT models, so skip these plots
ax = subplot(nPlotsPerModel,nModels,4);
title('Decisional choice model');
axis square;

ax = subplot(nPlotsPerModel,nModels,5);
title('Reaction time model');
axis square;

ax = subplot(nPlotsPerModel,nModels,6);
set(gca, 'FontName', 'Arial');
plot(ax,aloc,mean(featHemiPre),'LineWidth',0.25); hold on;
plot(ax,aloc,mean(featHemiPostVA),'LineWidth',0.25);
plot(ax,aloc,mean(featHemiPostAV),'LineWidth',0.25);
set(gca,'Xtick',aloc);
ylim([0.2,0.6]);
xlabel('Source location (degrees)');
ylabel(sprintf('First moment\n(mean activity)'));
axis square;

ax = subplot(nPlotsPerModel,nModels,7);
set(gca, 'FontName', 'Arial');
plot(ax,aloc,mean(featPlacePre),'LineWidth',0.25); hold on;
plot(ax,aloc,mean(featPlacePostVA),'LineWidth',0.25);
plot(ax,aloc,mean(featPlacePostAV),'LineWidth',0.25);
set(gca,'Xtick',aloc);
ylim([0,0.4]);
xlabel('Source location (degrees)');
axis square;

ax = subplot(nPlotsPerModel,nModels,8);
set(gca, 'FontName', 'Arial');
plot(ax,aloc,mean(featDecPre),'LineWidth',0.25); hold on;
plot(ax,aloc,mean(featDecPostVA),'LineWidth',0.25);
plot(ax,aloc,mean(featDecPostAV),'LineWidth',0.25);
set(gca,'Xtick',aloc);
ylim([0,1]);
xlabel('Source location (degrees)');
axis square

ax = subplot(nPlotsPerModel,nModels,9);
set(gca, 'FontName', 'Arial');
plot(ax,aloc,featPmfPre,'LineWidth',0.25); hold on;
plot(ax,aloc,featPmfPostVA,'LineWidth',0.25);
plot(ax,aloc,featPmfPostAV,'LineWidth',0.25);
set(gca,'Xtick',aloc);
ylim([0,1]);
ylabel('p responded right');
xlabel('Source location (degrees)');
axis square

ax = subplot(nPlotsPerModel,nModels,10);
set(gca, 'FontName', 'Arial');
plot(ax,aloc,featRtPre,'LineWidth',0.25); hold on;
plot(ax,aloc,featRtPostVA,'LineWidth',0.25);
plot(ax,aloc,featRtPostAV,'LineWidth',0.25);
set(gca,'Xtick',aloc);
% ylim([0,1]);
ylabel('Reaction time (s)');
xlabel('Source location (degrees)');
legend(conds);
axis square

ax = subplot(nPlotsPerModel,nModels,11);
set(gca, 'FontName', 'Arial');
RDM_hemi_recal = plotSecondMoment(ax,cat(2,featHemiPre,featHemiPostVA,featHemiPostAV)');
ylabel(sprintf('Second moment\n(pattern similarity)\nRecalibration'));

ax = subplot(nPlotsPerModel,nModels,12);
set(gca, 'FontName', 'Arial');
plotSecondMoment(ax,cat(2,featPlacePre,featPlacePostVA,featPlacePostAV)');

ax = subplot(nPlotsPerModel,nModels,13);
set(gca, 'FontName', 'Arial');
RDM_dec_recal = plotSecondMoment(ax,cat(2,featDecPre,featDecPostVA,featDecPostAV)');

ax = subplot(nPlotsPerModel,nModels,14);
set(gca, 'FontName', 'Arial');
RDM_pmf_recal = plotSecondMoment(ax,cat(2,featPmfPre,featPmfPostVA,featPmfPostAV)');

ax = subplot(nPlotsPerModel,nModels,15);
set(gca, 'FontName', 'Arial');
RDM_rt_recal = plotSecondMoment(ax,cat(2,featRtPre,featRtPostVA,featRtPostAV)');

ax = subplot(nPlotsPerModel,nModels,16);
set(gca, 'FontName', 'Arial');
plotSecondMoment(ax,repmat(featHemiPre,1,3)');
ylabel(sprintf('Second moment\n(pattern similarity)\nNo Recalibration'));

ax = subplot(nPlotsPerModel,nModels,17);
set(gca, 'FontName', 'Arial');
plotSecondMoment(ax,repmat(featPlacePre,1,3)');

ax = subplot(nPlotsPerModel,nModels,18);
set(gca, 'FontName', 'Arial');
plotSecondMoment(ax,repmat(featDecPre,1,3)');

ax = subplot(nPlotsPerModel,nModels,19);
set(gca, 'FontName', 'Arial');
plotSecondMoment(ax,repmat(featPmfPre,1,3)');

ax = subplot(nPlotsPerModel,nModels,20);
set(gca, 'FontName', 'Arial');
plotSecondMoment(ax,repmat(featRtPre,1,3)');

if plotCorr
    % Correlation plot of model predictions
    % First moment
    X = array2table(...
        [cat(2,mean(featHemiPre),mean(featHemiPostVA),mean(featHemiPostAV))',...
        cat(2,mean(featDecPre),mean(featDecPostVA),mean(featDecPostAV))',...
        cat(2,featPmfPre,featPmfPostVA,featPmfPostAV)',...
        cat(2,featRtPre,featRtPostVA,featRtPostAV)'],'VariableNames',...
        {'Sp','Dec','Pmf','Rt'});
    figure('Name','First moment correlation matrix');
    [R_1st,pval_1st] = corrplot(X,'type','Spearman','testR','on');
    % title('Correlation matrix of first moment predictions');
    % Second moment
    % Converting RDMs to lower triangular vectors
    temp = tril(RDM_hemi_recal);
    RDM_hemi_recal = temp(temp ~= 0);
    temp = tril(RDM_dec_recal);
    RDM_dec_recal = temp(temp ~= 0);
    temp = tril(RDM_pmf_recal);
    RDM_pmf_recal = temp(temp ~= 0);
    temp = tril(RDM_rt_recal);
    RDM_rt_recal = temp(temp ~= 0);
    X = array2table([tril(RDM_hemi_recal),RDM_dec_recal,RDM_pmf_recal,...
        RDM_rt_recal(:)],...
        'VariableNames',{'Sp','Dec','Pmf','Rt'});
    figure('Name','Second moment correlation matrix');
    [R_2nd,pval_2nd] = corrplot(X,'type','Spearman','testR','on');
end

% Writing tuning curve plot data to table
temp_model = repmat({'hemi'},numel(y_hemi),1);
temp_neuron = strcat('neuron_',cellfun(@num2str,(num2cell(1:size(y_hemi,1))),...
                                       'UniformOutput',false));
temp_neuron = repmat(temp_neuron,size(x_hemi,2),1);
temp_neuron = temp_neuron(:);
temp_aloc = repmat(x_hemi',size(y_hemi,1),1);
temp_act = y_hemi';
temp_act = temp_act(:);
tuning_hemi_table = table(temp_model,temp_neuron,temp_aloc,temp_act,'VariableNames',...
                          {'model','neuron','aloc','activation'});

temp_model = repmat({'place'},numel(y_place),1);
temp_neuron = strcat('neuron_',cellfun(@num2str,(num2cell(1:size(y_place,1))),...
                                       'UniformOutput',false));
temp_neuron = repmat(temp_neuron,size(x_place,2),1);
temp_neuron = temp_neuron(:);
temp_aloc = repmat(x_place',size(y_place,1),1);
temp_act = y_place';
temp_act = temp_act(:);
tuning_place_table = table(temp_model,temp_neuron,temp_aloc,temp_act,'VariableNames',...
                          {'model','neuron','aloc','activation'});

temp_model = repmat({'dec_unc'},numel(y_dec),1);
temp_neuron = strcat('neuron_',cellfun(@num2str,(num2cell(1:size(y_dec,1))),...
                                       'UniformOutput',false));
temp_neuron = repmat(temp_neuron,size(x_dec,2),1);
temp_neuron = temp_neuron(:);
temp_aloc = repmat(x_dec',size(y_dec,1),1);
temp_act = y_dec';
temp_act = temp_act(:);
tuning_dec_table = table(temp_model,temp_neuron,temp_aloc,temp_act,'VariableNames',...
                          {'model','neuron','aloc','activation'});
tuning_table = cat(1,tuning_hemi_table,tuning_place_table,tuning_dec_table);

% Writing mean tuning curves to table
conds = {'pre','postVA','postAV'};
models = {'hemi','place','dec_unc','dec_ch','rt'};
temp_aloc = repmat(aloc',numel(conds),1);
temp_aloc = repmat(temp_aloc,numel(models),1);
temp_conds = repmat(conds,numel(aloc),1);
temp_conds = temp_conds(:);
temp_conds = repmat(temp_conds,numel(models),1);
temp_model = repmat(models,numel(conds)*numel(aloc),1);
temp_model = temp_model(:);
temp_act_hemi = cat(2,mean(featHemiPre),mean(featHemiPostVA),mean(featHemiPostAV))';
temp_act_place = cat(2,mean(featPlacePre),mean(featPlacePostVA),mean(featPlacePostAV))';
temp_act_dec = cat(2,mean(featDecPre),mean(featDecPostVA),mean(featDecPostAV))';
temp_act_pmf = cat(2,featPmfPre,featPmfPostVA,featPmfPostAV)';
temp_act_rt = cat(2,featRtPre,featRtPostVA,featRtPostAV)';
mean_activation_table = table(temp_model,temp_conds,temp_aloc,...
                              cat(1,temp_act_hemi,temp_act_place,temp_act_dec,...
                                  temp_act_pmf,temp_act_rt),'VariableNames',...
                                  {'model','condition','aloc','mean_activation'});
varargout = {tuning_table,mean_activation_table};

if plotCorr
    varargout = cat(2,varargout,{R_1st,pval_1st,R_2nd,pval_2nd});
end

end

function [deg,y] = plotHemi(ax,propContra)

nNeurons = [25,25];
propContra = repmat(propContra,1,2);
popCenter = 90;
sigma = 67;
span = 20;
propIpsi = 1-propContra;
latRatio = 1;

deg = -180:180;

muContra = {};
muContra{1} = rand(floor(nNeurons(1)*propContra(1)),1)*span-(span/2)+popCenter;
muContra{2} = -(muContra{1}-popCenter)+popCenter-360;
muContra = cat(1,muContra{:});
mask = {};
mask{1} = false(floor(nNeurons(1)*propContra(1)),numel(deg));
mask{1}(:,ismember(deg,-180:-91)) = true;
mask{2} = false(floor(nNeurons(1)*propContra(1)),numel(deg));
mask{2}(:,ismember(deg,-89:180)) = true;

muIpsi = {};
muIpsi{1} = rand(ceil(nNeurons(1)*propIpsi(1)),1)*span-(span/2)-popCenter;
muIpsi{2} = -(muIpsi{1}+popCenter)-popCenter+360;
muIpsi = cat(1,muIpsi{:});
mask{3} = false(ceil(nNeurons(1)*propIpsi(1)),numel(deg));
mask{3}(:,ismember(deg,91:180)) = true;
mask{4} = false(ceil(nNeurons(1)*propIpsi(1)),numel(deg));
mask{4}(:,ismember(deg,-180:89)) = true;
mask = cat(1,mask{:});
                       
fun = @(m) normpdf(deg,m,sigma);

y = cat(1,arrayfun(fun,muContra,'UniformOutput',false),...
arrayfun(fun,muIpsi,'UniformOutput',false));
y = cat(1,y{:})/latRatio;
y = y./max(y(:));
y_plot = y;
y_plot(mask) = NaN;

plot(ax,deg,y_plot','Color',[0.5,0.5,0.5],'LineWidth',0.25);
set(gca,'Xtick',[-90,0,90]);
xlim([-180,180]);
axis square

end

function [deg,y] = plotDec(ax)
nNeurons = 5;
deg = -180:180;
% Decision model features
spanDec = 5;
sigmaDec = 10;
muDec = rand(nNeurons,1)*spanDec-(spanDec/2);
funDec = @(m,s) abs(abs(normcdf(deg,m+s,sigmaDec)-0.5)-0.5)*2;
y = arrayfun(funDec,muDec,zeros(size(muDec)),'UniformOutput',false);
y = cat(1,y{:});

plot(ax,deg,y','Color',[0.5,0.5,0.5],'LineWidth',0.25);
set(gca,'Xtick',[-90,0,90]);
xlim([-180,180]);
axis square

end

function [deg,y] = plotPlace(ax,mu)

sigma = 22.5;
latRatio = 1;

deg = -180:180;

% mu = (-270:15:270)';

fun = @(m) normpdf(deg,m,sigma);

y = arrayfun(fun,mu,'UniformOutput',false);
y = cat(1,y{:})/latRatio;
y = y./max(y(:));

plot(ax,deg,y','Color',[0.5,0.5,0.5],'LineWidth',0.25);
set(gca,'Xtick',[-90,0,90]);
xlim([-180,180]);
axis square
end

function RDM = plotSecondMoment(ax,Y)
plotMode = 'RDM';

G = Y*Y';
nConds = size(G,1);
H = eye(nConds)-ones(nConds)/nConds;
Gc = H*G*H';

if strcmp(plotMode,'G')
    imagesc(ax,Gc);
    set(gca,'XTick',[],'YTick',[]);
    axis square
    colormap('magma');
%     colorbar;
elseif strcmp(plotMode,'RDM')
    C = pcm_indicatorMatrix('allpairs',[1:nConds]');
    % Converting G to RDM
    RDM = squareform(diag(C'*Gc*C));
    RDM = rsa.util.rankTransform_equalsStayEqual(RDM);
    imagesc(ax,RDM);
    set(gca,'XTick',[],'YTick',[]);
    axis square
    colormap('magma');
%     colorbar;
end
end