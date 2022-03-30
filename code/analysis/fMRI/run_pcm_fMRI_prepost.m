function [out,M] = run_pcm_fMRI_prepost(varargin)

% Parsing input
p = inputParser;

validFitModes = {'group','individual'};
validFreeTypes = {'freechol','freedirect'};
validRunEffects = {'fixed','random'};
validPartVecModes = {'byRun','poolPrePost'};
validThirdModels = {'Pmf','Rt'};

addParameter(p,'fitRes',[],@(x) validateattributes(x,{'cell'},{'nonempty'}));
addParameter(p,'sigmaDec',10,@(x) validateattributes(x,{'numeric'},{'scalar','positive'}));
addParameter(p,'shift',2.3,@(x) validateattributes(x,{'numeric'},{'scalar','positive'}));
addParameter(p,'plotRes',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'plotIndiv',false,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'plotObs',false,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'centerG',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'topNvoxelTval',[],@(x) validateattributes(x,{'numeric'},{'scalar','positive'}));
addParameter(p,'topNvoxelSVM',[],@(x) validateattributes(x,{'numeric'},{'scalar','positive'}));
addParameter(p,'iseucnorm',false,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'fitMode','group',@(x) ismember(x,validFitModes));
addParameter(p,'runEffect','fixed',@(x) ismember(x,validRunEffects));
addParameter(p,'freeType','freechol',@(x) ismember(x,validFreeTypes));
addParameter(p,'partVecMode','byRun',@(x) ismember(x,validPartVecModes));
addParameter(p,'thirdModel','',@(x) ismember(x,validThirdModels));

parse(p,varargin{:});

fitRes = p.Results.fitRes; % Fitting results in case they are already computed
sigmaDec = p.Results.sigmaDec; % Sigma parameter of decision model
shift = p.Results.shift; % Shift parameter of all models, default derived from behavioural data
plotRes = p.Results.plotRes; % Weather to plot fitting results
plotIndiv = p.Results.plotIndiv; % Plot indiviudal fit results
plotObs = p.Results.plotObs; % Only plot the observed G matrixes
centerG = p.Results.centerG; % Weather to center model G matrix before fitting
topNvoxelTval = p.Results.topNvoxelTval; % Select top N voxels based on T values
topNvoxelSVM = p.Results.topNvoxelSVM; % Select top N voxels based on SVM
iseucnorm = p.Results.iseucnorm; % Weather to euclidean normalize observed data
fitMode = p.Results.fitMode; % Fitting mode
runEffect = p.Results.runEffect;
freeType = p.Results.freeType;
partVecMode = p.Results.partVecMode;
thirdModel = p.Results.thirdModel;

if plotObs && ~isempty(fitRes)
    error('Either the observed G matrices are plotted or the fit results, not both!');
elseif plotObs
    plotRes = false;
end

if strcmp(runEffect,'random') && centerG
    warning(['With runEffects set to ''random'' centering should not be done ',...
            'on the model G matrices. Setting centerG to false.']);
    centerG = false;
end

%% Getting fMRI data for pcm (nSubj x nROI)
if ~isempty(topNvoxelTval)
    [pcm_data,LUT] = pcm_getData(...
        'blocktype',{{'pretest','posttest-ladapt' 'posttest-radapt'}},...
        'output','bySession','topNvoxelTval',topNvoxelTval,'iseucnorm',iseucnorm);
elseif ~isempty(topNvoxelSVM)
    [pcm_data,LUT] = pcm_getData(...
        'blocktype',{{'pretest','posttest-ladapt' 'posttest-radapt'}},...
        'output','bySession','topNvoxelSVM',topNvoxelSVM,'iseucnorm',iseucnorm);
else
    [pcm_data,LUT] = pcm_getData(...
        'blocktype',{{'pretest','posttest-ladapt' 'posttest-radapt'}},...
        'output','bySession','iseucnorm',iseucnorm);
end

[nSubj,nROI] = size(pcm_data);
roiNames = {'HG','hA','IPL','IPS','FEF'};

%% Model specificaion
if ~ismember('fitRes',p.UsingDefaults)
    M = fitRes{2};
    out = fitRes{1};
else
    % Generating features
    aloc = [-12,-5,-2,0,2,5,12];
    nConds = numel(aloc)*3;
    
    % Hemifield model features
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
    fLoaded = load('fMRI_behav_all.mat');
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

    % Reaction time model features
    % To be added
    
    iModel = 1;
    % Order of 21 conditions: Pre,VA-post, AV-post
    % Model 1: Model with independent patterns (zero
    % correlation)
    M{iModel}.type       = 'component';
    M{iModel}.numGparams = 1;
    M{iModel}.Gc         = eye(nConds);
    M{iModel}.name       = 'null';
    iModel = iModel+1;
    
    % Hemifield with recalibration
    % Right hemisphere
    % Predicted response for each location for right hemifield neuron
    [Ac,Ac_pre,Ac_postVA,Ac_postAV] = deal([]);
    for i = 1:size(featHemiPre,1)
        % Predicted response for each location for right hemifield neurons
        Ac_pre(:,i,i)     = [featHemiPre(i,:),zeros(1,14)]';             % Pretest
        Ac_postVA(:,i,i)  = [zeros(1,7),featHemiPostVA(i,:),zeros(1,7)]'; % Recalibration VA
        Ac_postAV(:,i,i)  = [zeros(1,14),featHemiPostAV(i,:)]';           % Recalibration AV
    end
    Ac = cat(3,Ac_pre,Ac_postVA,Ac_postAV);
    theta0       = ones(size(Ac,3),1);
    Gc_hemi_recal      = getGmatrix(Ac,theta0,nConds,'center',centerG);
    
    % Hemifield no recalibration with two populations
    [Ac,Ac_pre,Ac_postVA,Ac_postAV] = deal([]);
    for i = 1:size(featHemiPre,1)
        % Predicted response for each location for right hemifield neurons
        Ac_pre(:,i,i)     = [featHemiPre(i,:),zeros(1,14)]';          % Pretest
        Ac_postVA(:,i,i)  = [zeros(1,7),featHemiPre(i,:),zeros(1,7)]'; % Recalibration VA No shift
        Ac_postAV(:,i,i)  = [zeros(1,14),featHemiPre(i,:)]';           % Recalibration AV No shift
    end
    Ac = cat(3,Ac_pre,Ac_postVA,Ac_postAV);
    theta0 = ones(size(Ac,3),1);
    Gc_hemi = getGmatrix(Ac,theta0,nConds,'center',centerG);
    
    % Decisional with recalibration
    [Ac,Ac_pre,Ac_postVA,Ac_postAV] = deal([]);
    for i = 1:size(featDecPre,1)
        % Predicted response for each location for decision model
        Ac_pre(:,i,i)     = [featDecPre(i,:),zeros(1,14)]';             % Pretest
        Ac_postVA(:,i,i)  = [zeros(1,7),featDecPostVA(i,:),zeros(1,7)]'; % Recalibration VA
        Ac_postAV(:,i,i)  = [zeros(1,14),featDecPostAV(i,:)]';           % Recalibration AV
    end
    Ac = cat(3,Ac_pre,Ac_postVA,Ac_postAV);
    theta0       = ones(size(Ac,3),1);
    Gc_dec_recal      = getGmatrix(Ac,theta0,nConds,'center',centerG);
    
    % Decisional no recalibration
    [Ac,Ac_pre,Ac_postVA,Ac_postAV] = deal([]);
    for i = 1:size(featDecPre,1)
        % Predicted response for each location for decision model
        Ac_pre(:,i,i)     = [featDecPre(i,:),zeros(1,14)]';           % Pretest
        Ac_postVA(:,i,i)  = [zeros(1,7),featDecPre(i,:),zeros(1,7)]'; % Recalibration VA NO shift
        Ac_postAV(:,i,i)  = [zeros(1,14),featDecPre(i,:)]';           % Recalibration AV No shift
    end
    Ac = cat(3,Ac_pre,Ac_postVA,Ac_postAV);
    theta0       = ones(size(Ac,3),1);
    Gc_dec      = getGmatrix(Ac,theta0,nConds,'center',centerG);
    
    if isempty(thirdModel)
        % Model 2: Hemifield no recal +Decisional no recal
        M{iModel}.type       = 'component';
        M{iModel}.numGparams = 2;
        M{iModel}.Gc(:,:,1)  = Gc_hemi;
        M{iModel}.Gc(:,:,2)  = Gc_dec;
        M{iModel}.name       = 'Hemi+Dec';
        iModel = iModel+1;
        
        % Model 3: Hemifield recal + Decisional
        M{iModel}.type       = 'component';
        M{iModel}.numGparams = 2;
        M{iModel}.Gc(:,:,1)  = Gc_hemi_recal;
        M{iModel}.Gc(:,:,2)  = Gc_dec;
        M{iModel}.name       = 'HemiRecal+Dec';
        iModel = iModel+1;
        
        % Model 4: Hemifield + Decisional recal
        M{iModel}.type       = 'component';
        M{iModel}.numGparams = 2;
        M{iModel}.Gc(:,:,1)  = Gc_hemi;
        M{iModel}.Gc(:,:,2)  = Gc_dec_recal;
        M{iModel}.name       = 'Hemi+DecRecal';
        iModel = iModel+1;
        
        % Model 5: Hemifield recal + Decisional recal
        M{iModel}.type       = 'component';
        M{iModel}.numGparams = 2;
        M{iModel}.Gc(:,:,1)  = Gc_hemi_recal;
        M{iModel}.Gc(:,:,2)  = Gc_dec_recal;
        M{iModel}.name       = 'HemiRecal+DecRecal';
        iModel = iModel+1;
        
    elseif strcmp(thirdModel,'Pmf')
        
        % Psychometric funtion with recalibration
        [Ac,Ac_pre,Ac_postVA,Ac_postAV] = deal([]);
        Ac_pre     = [featPmfPre,zeros(1,14)]';          % Pretest
        Ac_postVA  = [zeros(1,7),featPmfPostVA,zeros(1,7)]'; % Recalibration VA
        Ac_postAV  = [zeros(1,14),featPmfPostAV]';           % Recalibration AV
        Ac = cat(3,Ac_pre,Ac_postVA,Ac_postAV);
        theta0     = ones(size(Ac,3),1);
        Gc_pmf_recal     = getGmatrix(Ac,theta0,nConds,'center',centerG);
        
        % Psychometric funtion no recalibration
        [Ac,Ac_pre,Ac_postVA,Ac_postAV] = deal([]);
        Ac_pre     = [featPmfPre,zeros(1,14)]';          % Pretest
        Ac_postVA  = [zeros(1,7),featPmfPre,zeros(1,7)]'; % Recalibration VA
        Ac_postAV  = [zeros(1,14),featPmfPre]';           % Recalibration AV
        Ac = cat(3,Ac_pre,Ac_postVA,Ac_postAV);
        theta0     = ones(size(Ac,3),1);
        Gc_pmf = getGmatrix(Ac,theta0,nConds,'center',centerG);
        
        % Model 2: Hemifield + Dec + Pmf
        M{iModel}.type       = 'component';
        M{iModel}.numGparams = 3;
        M{iModel}.Gc(:,:,1)  = Gc_hemi;
        M{iModel}.Gc(:,:,2)  = Gc_dec;
        M{iModel}.Gc(:,:,3)  = Gc_pmf;
        M{iModel}.name       = 'Hemi+Dec+Pmf';
        iModel = iModel+1;
        
        % Model 3: Hemifield_recal + Dec + Pmf
        M{iModel}.type       = 'component';
        M{iModel}.numGparams = 3;
        M{iModel}.Gc(:,:,1)  = Gc_hemi_recal;
        M{iModel}.Gc(:,:,2)  = Gc_dec;
        M{iModel}.Gc(:,:,3)  = Gc_pmf;
        M{iModel}.name       = 'HemiRecal+Dec+Pmf';
        iModel = iModel+1;
        
        % Model 4: Hemifield + Dec_recal + Pmf
        M{iModel}.type       = 'component';
        M{iModel}.numGparams = 3;
        M{iModel}.Gc(:,:,1)  = Gc_hemi;
        M{iModel}.Gc(:,:,2)  = Gc_dec_recal;
        M{iModel}.Gc(:,:,3)  = Gc_pmf;
        M{iModel}.name       = 'Hemi+DecRecal+Pmf';
        iModel = iModel+1;
        
        % Model 5: Hemifield + Dec + Pmf_recal
        M{iModel}.type       = 'component';
        M{iModel}.numGparams = 3;
        M{iModel}.Gc(:,:,1)  = Gc_hemi;
        M{iModel}.Gc(:,:,2)  = Gc_dec;
        M{iModel}.Gc(:,:,3)  = Gc_pmf_recal;
        M{iModel}.name       = 'Hemi+Dec+PmfRecal';
        iModel = iModel+1;
        
        % Model 6: Hemifield + Dec_recal + Pmf_recal
        M{iModel}.type       = 'component';
        M{iModel}.numGparams = 3;
        M{iModel}.Gc(:,:,1)  = Gc_hemi;
        M{iModel}.Gc(:,:,2)  = Gc_dec_recal;
        M{iModel}.Gc(:,:,3)  = Gc_pmf_recal;
        M{iModel}.name       = 'Hemi+DecRecal+PmfRecal';
        iModel = iModel+1;
        
        % Model 7: Hemifield_recal + Dec + Pmf_recal
        M{iModel}.type       = 'component';
        M{iModel}.numGparams = 3;
        M{iModel}.Gc(:,:,1)  = Gc_hemi_recal;
        M{iModel}.Gc(:,:,2)  = Gc_dec;
        M{iModel}.Gc(:,:,3)  = Gc_pmf_recal;
        M{iModel}.name       = 'HemiRecal+Dec+PmfRecal';
        iModel = iModel+1;
        
        % Model 8: Hemifield_recal + Dec_recal + Pmf
        M{iModel}.type       = 'component';
        M{iModel}.numGparams = 3;
        M{iModel}.Gc(:,:,1)  = Gc_hemi_recal;
        M{iModel}.Gc(:,:,2)  = Gc_dec_recal;
        M{iModel}.Gc(:,:,3)  = Gc_pmf;
        M{iModel}.name       = 'HemiRecal+DecRecal+Pmf';
        iModel = iModel+1;
        
        % Model 9: Hemifield_recal + Dec_recal + Pmf_recal
        M{iModel}.type       = 'component';
        M{iModel}.numGparams = 3;
        M{iModel}.Gc(:,:,1)  = Gc_hemi_recal;
        M{iModel}.Gc(:,:,2)  = Gc_dec_recal;
        M{iModel}.Gc(:,:,3)  = Gc_pmf_recal;
        M{iModel}.name       = 'HemiRecal+DecRecal+PmfRecal';
        iModel = iModel+1;
        
    elseif strcmp(thirdModel,'Rt')
        error('PCM with reaction time model not implemented. ');
    end
    
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
        out = cell(nROI,1);
        parfor iROI = 1:nROI
            Y = pcm_data(:,iROI)';
            partVec = cellfun(@(x) getPartVec(x,partVecMode),LUT(:,iROI),...
                              'UniformOutput',false);
            % Converting aloc to integer conditions
            condVec = cellfun(@(x) x.aloc,LUT(:,iROI),'UniformOutput',false);
            isRight = cellfun(@(x) double(ismember(x.blocktype,'posttest-radapt')),...
                LUT(:,iROI),'UniformOutput',false);
            isLeft = cellfun(@(x) double(ismember(x.blocktype,'posttest-ladapt')),...
                LUT(:,iROI),'UniformOutput',false);
            [~,~,condVec] = cellfun(@unique,condVec,'UniformOutput',false);
            condVec = cellfun(@(x,y,z) x+(y*7)+(z*14),condVec,isLeft,isRight,'UniformOutput',false);
            if strcmp(fitMode,'group')
                y = cell(1,6);
            else
                y = cell(1,5);
            end
            [y{:}] = pcm_fit(Y,M,condVec,partVec,'fitMode',fitMode,...
                                'runEffect',runEffect);
            out{iROI} = y;
        end
        out = cat(1,out{:});
    else
        % Just plot the observed G matrices (as they go in the analysis)
        hFig = arrayfun(@(x) figure(),1:(size(pcm_data,1)+1));
        for iROI = 1:nROI
            Y = pcm_data(:,iROI)';
            partVec = cellfun(@(x) getPartVec(x,partVecMode),LUT(:,iROI),....
                              'UniformOutput',false);
            % Converting aloc to integer conditions
            condVec = cellfun(@(x) x.aloc,LUT(:,iROI),'UniformOutput',false);
            isRight = cellfun(@(x) double(ismember(x.blocktype,'posttest-radapt')),...
                LUT(:,iROI),'UniformOutput',false);
            isLeft = cellfun(@(x) double(ismember(x.blocktype,'posttest-ladapt')),...
                LUT(:,iROI),'UniformOutput',false);
            [~,~,condVec] = cellfun(@unique,condVec,'UniformOutput',false);
            condVec = cellfun(@(x,y,z) x+(y*7)+(z*14),condVec,isLeft,isRight,'UniformOutput',false);
            plotObservedG(Y,condVec,partVec,runEffect,iROI,roiNames{iROI},hFig);
        end
    end
end

%% Plot fit results
if plotRes
    % Plot all models except baseline and free
    modelsToPlot = 2:(iModel-1);
    if plotIndiv
        for iSub = 1:size(pcm_data,1)
            hFig = [figure(),figure(),figure()];
            for iROI = 1:nROI
                Y = pcm_data(:,iROI)';
                partVec = cellfun(@(x) getPartVec(x,partVecMode),LUT(:,iROI),...
                                  'UniformOutput',false);
                % Converting aloc to integer conditions
                condVec = cellfun(@(x) x.aloc,LUT(:,iROI),'UniformOutput',false);
                isRight = cellfun(@(x) double(ismember(x.blocktype,'posttest-radapt')),...
                    LUT(:,iROI),'UniformOutput',false);
                isLeft = cellfun(@(x) double(ismember(x.blocktype,'posttest-ladapt')),...
                    LUT(:,iROI),'UniformOutput',false);
                [~,~,condVec] = cellfun(@unique,condVec,'UniformOutput',false);
                condVec = cellfun(@(x,y,z) x+(y*7)+(z*14),condVec,isLeft,isRight,'UniformOutput',false);
                pcm_plot_fitRes(out,M,Y,condVec,partVec,runEffect,iROI,modelsToPlot,...
                    hFig,roiNames,'iSub',iSub);
            end
            
            % Setting the scale of all theta plots' y axes uniform
            ax = findobj(hFig(1),'Type','Axes');
            y = get(ax(1:3:15),'ylim');
            y = cat(1,y{:});
            yl = [min(y(:)),max(y(:))];
            linkaxes(ax(1:3:15),'y');
            ylim(ax(1),yl);
            
            % Setting super title
            if iROI == nROI
                for iFig = 1:3
                    set(0, 'currentfigure', hFig(iFig));
                    suplabel(sprintf('sub-%d',iSub),'t',[.08 .08 .84 .86]);
                end
            end
        end
    else
        hFig = [figure(),figure(),figure()];
        for iROI = 1:nROI
            Y = pcm_data(:,iROI)';
            partVec = cellfun(@(x) getPartVec(x,partVecMode),LUT(:,iROI),...
                              'UniformOutput',false);
            % Converting aloc to integer conditions
            condVec = cellfun(@(x) x.aloc,LUT(:,iROI),'UniformOutput',false);
            isRight = cellfun(@(x) double(ismember(x.blocktype,'posttest-radapt')),...
                LUT(:,iROI),'UniformOutput',false);
            isLeft = cellfun(@(x) double(ismember(x.blocktype,'posttest-ladapt')),...
                LUT(:,iROI),'UniformOutput',false);
            [~,~,condVec] = cellfun(@unique,condVec,'UniformOutput',false);
            condVec = cellfun(@(x,y,z) x+(y*7)+(z*14),condVec,isLeft,isRight,'UniformOutput',false);
            pcm_plot_fitRes(out,M,Y,condVec,partVec,runEffect,iROI,modelsToPlot,hFig,roiNames);
        end
        
        % Setting the scale of all theta plots' y axes uniform
        ax = findobj(hFig(1),'Type','Axes');
        y = get(ax(1:3:15),'ylim');
        y = cat(1,y{:});
        yl = [min(y(:)),max(y(:))];
        linkaxes(ax(1:3:15),'y');
        ylim(ax(1),yl);
    end
end

end

function partVec = getPartVec(LUT,partVecMode)

switch partVecMode
    case 'poolPrePost'
        un = unique(LUT(:,{'session','blocktype'}),'rows');
        idxPostL = un.session(ismember(un.blocktype,'posttest-ladapt'));
        idxPostR = un.session(ismember(un.blocktype,'posttest-radapt'));
        partVec = NaN(size(LUT,1),1);
        for i = 1:numel(idxPostL)
            partVec(ismember(LUT.session,[idxPostL(i),idxPostR(i)])) = i;
        end
    case 'byRun'
        partVec = LUT.spm_session;
    otherwise
        error('Unrecognized value for partVecMode!');
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
    ax = subplot(1,5,iROI);
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
    if iROI == 5
        suplabel(sprintf('sub-%d',s),'t',[.08 .08 .84 .86]);
    end
end
% Plot subject level G matrices
set(0,'currentfigure',hFig(s+1));
set(gcf,'Position',[300 40 1400 320]);
ax = subplot(1,5,iROI);
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
if iROI == 5
    suplabel('group','t',[.08 .08 .84 .86]);
end
end