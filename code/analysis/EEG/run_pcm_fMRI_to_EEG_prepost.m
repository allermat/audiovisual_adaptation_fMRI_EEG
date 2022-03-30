function [out,M] = run_pcm_fMRI_to_EEG_prepost(varargin)

% Parsing input
p = inputParser;

validRunEffects = {'fixed','random'};
validFitModes = {'group','individual'};
validErpModes = {'byRun','rand','glmByRun'};
validFreeTypes = {'freechol','freedirect'};

addParameter(p,'fitRes',[],@(x) validateattributes(x,{'cell'},{'nonempty'}));
addParameter(p,'plotRes',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'centerG',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'runEffect','fixed',@(x) ismember(x,validRunEffects));
addParameter(p,'fitMode','group',@(x) ismember(x,validFitModes));
addParameter(p,'erpMode','byRun',@(x) ismember(x,validErpModes));
addParameter(p,'freeType','freechol',@(x) ismember(x,validFreeTypes));
addParameter(p,'iseucnorm',false,@(x) validateattributes(x,{'logical'},{'scalar'}));

parse(p,varargin{:});

fitRes = p.Results.fitRes;
plotRes = p.Results.plotRes;
centerG = p.Results.centerG; % Weather to center model G matrix before fitting
runEffect = p.Results.runEffect;
fitMode = p.Results.fitMode; % Fitting mode
erpMode = p.Results.erpMode; % How the compute ERPs for the analysis
freeType = p.Results.freeType;
iseucnorm = p.Results.iseucnorm; % Weather to euclidean normalize observed data

if strcmp(runEffect,'random') && centerG
    warning(['With runEffects set to ''random'' centering should not be done ',...
            'on the model G matrices. Setting centerG to false.']);
    centerG = false;
end

%% Getting data for pcm 

% Getting data from EEG (nSubj x nTimePoints)
timeWins = cellfun(@(x) (0:5:100)+x,num2cell(50:100:350),'UniformOutput',false);
% timeWins = {50:5:150,150:5:450};
% timeWins = cellfun(@(x) (0:5:50)+x,num2cell(50:50:400),'UniformOutput',false);
% timeWins = cellfun(@(x) (0:5:25)+x,num2cell(50:25:425),'UniformOutput',false);
temp = cellfun(@(x) num2cell(x([1,end])'),timeWins,'UniformOutput',false);
temp = cellfun(@num2str,cat(2,temp{:}),'UniformOutput',false);
timeWinStr = strcat(temp(1,:),'-',temp(2,:),' ms');
[pcm_data_EEG,LUT_EEG] = pcm_getData_EEG(timeWins,erpMode,'iseucnorm',iseucnorm);
nConds = numel(unique(LUT_EEG{1}.rdmGrVar));
nTimeWin = numel(timeWinStr);
    
%% Model specificaion
if ~ismember('fitRes',p.UsingDefaults)
    M = fitRes{2};
    out = fitRes{1};
else
    % Using fMRI data as components to explain EEG activation patterns
    % fMRI (nSubj x nROI)
    % IF 'byHemisphere',true is specified, then each ROI is further
    % subdivided to two hemishperes, 
    % order: ROI1_lh, ROI2_lh...,ROI1_rh...ROIn_rh
    [pcm_data_fMRI,LUT_fMRI] = pcm_getData(...
        'blocktype',{{'pretest','posttest-ladapt' 'posttest-radapt'}},...
        'output','bySession','byHemisphere',true,'iseucnorm',iseucnorm);
    % Swapping IPS and IPL
    pcm_data_fMRI = pcm_data_fMRI(:,[1,2,4,3,5:7,9,8,10]);
    [nSubj,nROI] = size(pcm_data_fMRI);
    roiNames = {'HG','hA','IPS','IPL','FEF'};
    Gc_fMRI = cell(nROI,1);
    for iROI = 1:nROI
        Y = pcm_data_fMRI(:,iROI)';
        partVec = cellfun(@(x) x.spm_session,LUT_fMRI(:,iROI),'UniformOutput',false);
        % Converting aloc to integer conditions
        condVec = cellfun(@(x) x.aloc,LUT_fMRI(:,iROI),'UniformOutput',false);
        isRight = cellfun(@(x) double(ismember(x.blocktype,'posttest-radapt')),...
            LUT_fMRI(:,iROI),'UniformOutput',false);
        isLeft = cellfun(@(x) double(ismember(x.blocktype,'posttest-ladapt')),...
            LUT_fMRI(:,iROI),'UniformOutput',false);
        [~,~,condVec] = cellfun(@unique,condVec,'UniformOutput',false);
        condVec = cellfun(@(x,y,z) x+(y*7)+(z*14),condVec,isLeft,isRight,'UniformOutput',false);
        Gc_fMRI{iROI} = getGmatrix_fMRI(Y,partVec,condVec,centerG);
    end
    
    % Order of 21 conditions: Pre,VA-post, AV-post
    % Model 1: Model with independent patterns (zero
    % correlation)
    M{1}.type       = 'component';
    M{1}.numGparams = 1;
    M{1}.Gc         = eye(nConds);
    M{1}.name       = 'null';
    
    % Models for the ROIs
    for i = 1:numel(roiNames)
        M{i+1}.type       = 'component';
        M{i+1}.numGparams = 2;
        M{i+1}.Gc(:,:,1)  = Gc_fMRI{i};
        M{i+1}.Gc(:,:,2)  = Gc_fMRI{i+5};
        M{i+1}.name       = roiNames{i};
    end
    
%     % ROIs combined
%     M{i+2}.type       = 'component';
%     M{i+2}.numGparams = 4;
%     M{i+2}.Gc(:,:,1)  = Gc_fMRI{1};
%     M{i+2}.Gc(:,:,2)  = Gc_fMRI{6};
%     M{i+2}.Gc(:,:,3)  = Gc_fMRI{2};
%     M{i+2}.Gc(:,:,4)  = Gc_fMRI{7};
%     M{i+2}.name       = [roiNames{1},'-',roiNames{2}];
%     
%     M{i+3}.type       = 'component';
%     M{i+3}.numGparams = 4;
%     M{i+3}.Gc(:,:,1)  = Gc_fMRI{3};
%     M{i+3}.Gc(:,:,2)  = Gc_fMRI{8};
%     M{i+3}.Gc(:,:,3)  = Gc_fMRI{4};
%     M{i+3}.Gc(:,:,4)  = Gc_fMRI{9};
%     M{i+3}.name       = [roiNames{3},'-',roiNames{4}];
%     
%     M{i+4}.type       = 'component';
%     M{i+4}.numGparams = 6;
%     M{i+4}.Gc(:,:,1)  = Gc_fMRI{3};
%     M{i+4}.Gc(:,:,2)  = Gc_fMRI{8};
%     M{i+4}.Gc(:,:,3)  = Gc_fMRI{4};
%     M{i+4}.Gc(:,:,4)  = Gc_fMRI{9};
%     M{i+4}.Gc(:,:,5)  = Gc_fMRI{5};
%     M{i+4}.Gc(:,:,6)  = Gc_fMRI{10};
%     M{i+4}.name       = [roiNames{3},'-',roiNames{4},'-',roiNames{5}];
    
%     M{i+2}.type       = 'component';
%     M{i+2}.numGparams = nROI;
%     for j = 1:nROI
%         M{i+2}.Gc(:,:,j)  = Gc_fMRI{j};
%     end
%     M{i+2}.name       = 'All';
    
    % Free model as Noise ceiling
    if strcmp(freeType,'freechol')
        M{i+2}.type       = 'freechol';
        M{i+2}.numCond    = nConds;
        M{i+2}.name       = 'noiseceiling';
        M{i+2}            = pcm_prepFreeModel(M{i+2});
    else
        M{i+2}.type       = 'freedirect';
        M{i+2}.numGparams = 0;
        M{i+2}.name       = 'noiseceiling';
        M{i+2}.theta0     = [];
    end
    
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
        pcm_plot_fitRes(out,M,Y,condVec,partVec,runEffect,iTw,modelsToPlot,hFig,timeWinStr,...
            'relTo',relTo);
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