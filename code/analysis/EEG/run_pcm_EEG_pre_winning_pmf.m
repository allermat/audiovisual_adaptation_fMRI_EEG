function [out,M] = run_pcm_EEG_pre_winning_pmf(varargin)

% Parsing input
p = inputParser;

validRunEffects = {'fixed','random'};
validFitModes = {'group','individual'};
validErpModes = {'byRun','rand','glmByRun'};
validFreeTypes = {'freechol','freedirect'};

addParameter(p,'fitRes',[],@(x) validateattributes(x,{'cell'},{'nonempty'}));
addParameter(p,'plotRes',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'plotObs',false,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'centerG',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'runEffect','fixed',@(x) ismember(x,validRunEffects));
addParameter(p,'fitMode','group',@(x) ismember(x,validFitModes));
addParameter(p,'erpMode','byRun',@(x) ismember(x,validErpModes));
addParameter(p,'freeType','freechol',@(x) ismember(x,validFreeTypes));
addParameter(p,'sigmaDec',10,@(x) validateattributes(x,{'numeric'},{'scalar','positive'}));
addParameter(p,'iseucnorm',false,@(x) validateattributes(x,{'logical'},{'scalar'}));

parse(p,varargin{:});

fitRes = p.Results.fitRes;
plotRes = p.Results.plotRes;
plotObs = p.Results.plotObs;
centerG = p.Results.centerG; % Weather to center model G matrix before fitting
runEffect = p.Results.runEffect;
fitMode = p.Results.fitMode; % Fitting mode
erpMode = p.Results.erpMode; % How the compute ERPs for the analysis
freeType = p.Results.freeType;
sigmaDec = p.Results.sigmaDec;
iseucnorm = p.Results.iseucnorm;

if strcmp(runEffect,'random') && centerG
    warning(['With runEffects set to ''random'' centering should not be done ',...
            'on the model G matrices. Setting centerG to false.']);
    centerG = false;
end

if plotObs && ~isempty(fitRes)
    error('Either the observed G matrices are plotted or the fit results, not both!');
elseif plotObs
    plotRes = false;
end

%% Getting data for pcm 

% Getting data from EEG (nSubj x nTimePoints)
timeWins = cellfun(@(x) (0:5:100)+x,num2cell(50:100:350),'UniformOutput',false);
temp = cellfun(@(x) num2cell(x([1,end])'),timeWins,'UniformOutput',false);
temp = cellfun(@num2str,cat(2,temp{:}),'UniformOutput',false);
timeWinStr = strcat(temp(1,:),'-',temp(2,:),' ms');
[pcm_data_EEG,LUT_EEG] = pcm_getData_EEG(timeWins,erpMode,'blocktype','pre',...
                                         'iseucnorm',iseucnorm);
nConds = numel(unique(LUT_EEG{1}.rdmGrVar));
nTimeWin = numel(timeWinStr);
    
%% Model specificaion
if ~ismember('fitRes',p.UsingDefaults)
    M = fitRes{2};
    out = fitRes{1};
else   
    % Generating features for different models
    aloc = [-12,-5,-2,0,2,5,12];
    shift = 2.3; % Derived from behavioural data
    sigmaHemi = 64; % Based on Salminen et al. 2009
    popCenter = 90;
    nNeurons = 360;
    propR = 0.7;
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

    % Decisional uncertainty model features
    spanDec = 2;
    muDec = rand(nNeurons,1)*spanDec-(spanDec/2);
    funDec = @(m,s) abs(abs(normcdf(aloc,m+s,sigmaDec)-0.5)-0.5)*2;
    y = arrayfun(funDec,muDec,zeros(size(muDec)),'UniformOutput',false);
    featDecPre = cat(1,y{:});
    y = arrayfun(funDec,muDec,shift*ones(size(muDec)),'UniformOutput',false);
    featDecPostVA = cat(1,y{:});
    y = arrayfun(funDec,muDec,-shift*ones(size(muDec)),'UniformOutput',false);
    featDecPostAV = cat(1,y{:});
    
    % Psychometric funtion model features
    % Load PMF data
    fLoaded = load('EEG_behav_all.mat');
    pre = cat(2, fLoaded.out.pre);
    post = cat(2, fLoaded.out.post);
    % Computing group average PMF values
    y = arrayfun(@(i) mean(pre(i).NumPos./pre(i).OutOfNum), 1:size(pre,2),...
                 'UniformOutput', false);
    featPmfPre = mean(cat(1,y{:}));
    % post-adaptation
    [y_AV,y_VA] = deal({});
    for i = 1:size(post,2)
        idx_AV = 1:(size(post(i).NumPos,1)/2);
        idx_VA = ((size(post(i).NumPos,1)/2)+1):size(post(i).NumPos,1);
        y_AV{i} = mean(post(i).NumPos(idx_AV,:)./post(i).OutOfNum(idx_AV,:),1);
        y_VA{i} = mean(post(i).NumPos(idx_VA,:)./post(i).OutOfNum(idx_VA,:),1);
    end
    % post-AV adaptation
    featPmfPostAV = mean(cat(1,y_AV{:}));
    % post-VA adaptation
    featPmfPostVA = mean(cat(1,y_VA{:}));

    % Reaction time model features
    % To be added
    
    iModel = 1;
    % Model 1: Model with independent patterns (zero
    % correlation)
    M{iModel}.type       = 'component';
    M{iModel}.numGparams = 1;
    M{iModel}.Gc         = eye(nConds);
    M{iModel}.name       = 'null';
    iModel = iModel+1;
    
    % Model 2: Hemifield no recalibration
    % Right hemisphere
    % Predicted response for each location for right hemifield neuron
    Ac = [];
    for i = 1:size(featHemiPre,1)
        % Predicted response for each location for right hemifield neurons
        Ac(:,i,i)     = featHemiPre(i,:)';             % Pretest
    end
    theta0 = ones(size(Ac,3),1);
    Gc_hemi = getGmatrix(Ac,theta0,nConds,'center',centerG);
    
    M{iModel}.type       = 'component';
    M{iModel}.numGparams = 1;
    M{iModel}.Gc         = Gc_hemi;
    M{iModel}.name       = 'Hemi';
    iModel = iModel+1;
    
    % Decisional no recalibration
    Ac = [];
    for i = 1:size(featDecPre,1)
        % Predicted response for each location for decision model
        Ac(:,i,i) = featDecPre(i,:)';             % Pretest
    end
    theta0 = ones(size(Ac,3),1);
    Gc_dec = getGmatrix(Ac,theta0,nConds,'center',centerG);
    
    % Psychometric funtion no recalibration
    Ac = featPmfPre';          % Pretest
    theta0     = ones(size(Ac,3),1);
    Gc_pmf = getGmatrix(Ac,theta0,nConds,'center',centerG);
    
    % Model 3: Hemifield+Decisional
    M{iModel}.type       = 'component';
    M{iModel}.numGparams = 2;
    M{iModel}.Gc(:,:,1)  = Gc_hemi;
    M{iModel}.Gc(:,:,2)  = Gc_dec;
    M{iModel}.name       = 'Hemi+Dec';
    iModel = iModel+1;
    
    % Model 4: Hemifield+Pmf
    M{iModel}.type       = 'component';
    M{iModel}.numGparams = 2;
    M{iModel}.Gc(:,:,1)  = Gc_hemi;
    M{iModel}.Gc(:,:,2)  = Gc_pmf;
    M{iModel}.name       = 'Hemi+Pmf';
    iModel = iModel+1;
    
    % Model 5: Hemifield+Decisional+Pmf
    M{iModel}.type       = 'component';
    M{iModel}.numGparams = 3;
    M{iModel}.Gc(:,:,1)  = Gc_hemi;
    M{iModel}.Gc(:,:,2)  = Gc_dec;
    M{iModel}.Gc(:,:,3)  = Gc_pmf;
    M{iModel}.name       = 'Hemi+Dec+Pmf';
    iModel = iModel+1;

    % Last model: Free model as Noise ceiling
    if strcmp(freeType,'freechol')
        M{iModel}.type      = 'freechol';
        M{iModel}.numCond   = nConds;
        M{iModel}.name      = 'noiseceiling';
        M{iModel}           = pcm_prepFreeModel(M{iModel});
    else
        M{iModel}.type       = 'freedirect';
        M{iModel}.numGparams = 0;
        M{iModel}.name       = 'noiseceiling';
        M{iModel}.theta0     = [];
    end
    
    if ~plotObs
        % Fitting models
        out = cell(nTimeWin,1);
        parfor iTw = 1:nTimeWin
            Y = pcm_data_EEG(:,iTw)';
            partVec = cellfun(@(x) x.partVec,LUT_EEG(:,iTw),'UniformOutput',false);
            condVec = cellfun(@(x) x.rdmGrVar,LUT_EEG(:,iTw),'UniformOutput',false);
            if strcmp(fitMode,'group')
                temp = cell(1,6);
            else
                temp = cell(1,5);
            end
            [temp{:}] = pcm_fit(Y,M,condVec,partVec,'fitMode',fitMode,...
                'runEffect',runEffect);
            out{iTw} = temp;
        end
        out = cat(1,out{:});
    else
        % Just plot the observed G matrices (as they go in the analysis)
        hFig = arrayfun(@(x) figure(),1:(size(pcm_data_EEG,1)+1));
        for iTw = 1:nTimeWin
            Y = pcm_data_EEG(:,iTw)';
            partVec = cellfun(@(x) x.partVec,LUT_EEG(:,iTw),'UniformOutput',false);
            % Converting aloc to integer conditions
            condVec = cellfun(@(x) x.rdmGrVar,LUT_EEG(:,iTw),'UniformOutput',false);
            plotObservedG(Y,condVec,partVec,runEffect,iTw,timeWinStr{iTw},hFig);
        end
    end
end

%% Plot fit results
if plotRes
    hFig = [figure(),figure(),figure()];
    % Plotting relative to this model
    relTo = 2;
    modelsToPlot = 2:(numel(M)-1);
    for iTw = 1:nTimeWin
        Y = pcm_data_EEG(:,iTw)';
        partVec = cellfun(@(x) x.partVec,LUT_EEG(:,iTw),'UniformOutput',false);
        condVec = cellfun(@(x) x.rdmGrVar,LUT_EEG(:,iTw),'UniformOutput',false);
        pcm_plot_fitRes(out,M,Y,condVec,partVec,runEffect,iTw,modelsToPlot,hFig,...
                        timeWinStr,'relTo',relTo);
    end
    
    % Setting the scale of all theta plots' y axes uniform
    ax = findobj(hFig(1),'Type','Axes');
    temp = get(ax(1:3:(nTimeWin*3)),'ylim');
    temp = cat(1,temp{:});
    yl = [min(temp(:)),max(temp(:))];
    linkaxes(ax(1:3:(nTimeWin*3)),'y');
    ylim(ax(1),yl);
end

end

function plotObservedG(Y,conditionVec,partitionVec,runEffect,iROI,roiName,hFig)
% Plot observed G matrices as they go in the PCM analysis. 
% This is code from pcm_setUpFit.m where the crossvalidated G matrix is
% estimated between lines 56-98. 

doCenter = false;
nSubj = numel(Y);
nConds = numel(unique(conditionVec{1}));
G_hat = NaN(nConds,nConds,nSubj);
H = eye(nConds)-ones(nConds)/nConds;
for s = 1:nSubj
    % Set up the main matrices
    [N(s,1),P(s,1)] = size(Y{s});  
    
    cV = conditionVec{s};
    pV = partitionVec{s};
    
    switch runEffect
        case 'random'
            B{s}   = pcm_indicatorMatrix('identity_p',pV);
            X{s}   = zeros(N(s),0);
            numPart=size(B{s},2);
            run0(s,1)=log(sum(sum((pinv(B{s})*Y{s}).^2))./(numPart*P(s)));
            RX = eye(N(s))-B{s}*pinv(B{s});
            G_hat(:,:,s) = pcm_estGCrossval(Y{s},pV,cV);
        case 'fixed'
            B{s}  =  zeros(N(s),0);
            run0  =  [];
            X{s}  =  pcm_indicatorMatrix('identity_p',pV);
            RX = eye(N(s))-X{s}*pinv(X{s});
            G_hat(:,:,s) = pcm_estGCrossval(RX*Y{s},pV,cV);
    end
    % Center G matrices if necessary
    if doCenter
        G_toPlot = H*G_hat(:,:,s)*H';
    else
        G_toPlot = G_hat(:,:,s);
    end
    % Plot subject level G matrices
    set(0,'currentfigure',hFig(s));
    set(gcf,'Position',[300 40 1400 320]);
    ax = subplot(1,4,iROI);
    imagesc(ax,G_toPlot);
    set(gca,'XTick',[],'YTick',[]);
    axis square
    colorbar;
%     colormap('magma')
    title(roiName);
    % If there are multiples of the 7 spatial locations in the covariance
    % matrix, then plot a grid on top for easier orientation. 
    grid = floor(size(G_hat,1)/7); 
    if grid > 1
        % Computing grid line coordinates
        width = size(G_hat,1);
        x = (0:round(width/grid):width)+0.5;
        x([1,end]) = [];
        x = repmat(x,2,1);
        y = repmat([0,width+0.5]',1,2);
        % Drawing vertical grid lines
        line(x,y,'Color','w','LineWidth',1);
        % Drawing horizontal grid lines
        line(y,x,'Color','w','LineWidth',1);
    end
    if iROI == 4
        suplabel(sprintf('sub-%d',s),'t',[.08 .08 .84 .86]);
    end
end
% Plot subject level G matrices
set(0,'currentfigure',hFig(s+1));
set(gcf,'Position',[300 40 1400 320]);
ax = subplot(1,4,iROI);
if doCenter
    imagesc(ax,H*mean(G_hat,3)*H');
else
    imagesc(ax,mean(G_hat,3));
end
set(gca,'XTick',[],'YTick',[]);
axis square
colorbar
% colormap('magma')
title(roiName);
% If there are multiples of the 7 spatial locations in the covariance
% matrix, then plot a grid on top for easier orientation. 
if grid > 1
    % Computing grid line coordinates
    width = size(G_hat,1);
    x = (0:round(width/grid):width)+0.5;
    x([1,end]) = [];
    x = repmat(x,2,1);
    y = repmat([0,width+0.5]',1,2);
    % Drawing vertical grid lines
    line(x,y,'Color','w','LineWidth',1);
    % Drawing horizontal grid lines
    line(y,x,'Color','w','LineWidth',1);
end
% Title above all subplots
if iROI == 4
    suplabel('group','t',[.08 .08 .84 .86]);
end
end

function Gc = getGmatrix_fMRI(Y,partVec,condVec,doCenter)

for s=1:length(Y)
    G_hat(:,:,s)=pcm_estGCrossval(Y{s},partVec{s},condVec{s}); 
end
Gm = mean(G_hat,3); % Mean estimate  
nConds = numel(unique(condVec{s}));
if doCenter
    H = eye(nConds)-ones(nConds)/nConds; 
    Gc = H*Gm*H';
else
    Gc = Gm;
end

end

function Gc = getGmatrix_EEG(Y)

for s=1:length(Y)
    G_hat(:,:,s) = Y{s}*Y{s}';
end
G = mean(G_hat,3);
nConds = size(G,1);
H = eye(nConds)-ones(nConds)/nConds;
Gc = H*G*H';

end