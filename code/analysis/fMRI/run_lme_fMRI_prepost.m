function varargout = run_lme_fMRI_prepost(varargin)

% Parsing input
p = inputParser;

validThirdModels = {'Pmf','Rt'};
validHemispheres = {'lh','rh','mean_lhrh','lh_and_rh'};

addParameter(p,'hemisphere','lh_and_rh',@(x) ismember(x,validHemispheres));
addParameter(p,'thirdModel','',@(x) ismember(x,validThirdModels));
addParameter(p,'plotFigure',false,@(x) validateattributes(x,{'logical'},{'scalar'}));

parse(p,varargin{:});

hemisphere = p.Results.hemisphere;
thirdModel = p.Results.thirdModel;
plotFigure = p.Results.plotFigure;

plotModelDesign = false;

if strcmp(hemisphere,'lh_and_rh')
    hemiToLoad = {'lh','rh'};
else
    hemiToLoad = {hemisphere};
end
if isempty(thirdModel)
    modelsToUse = {'Const','Sp','Dec'};
elseif strcmp(thirdModel,'Pmf')
    modelsToUse = {'Const','Sp','Dec','Pmf'};
else
    error('Not yet implemented!');
end
[BIC,AIC,LL] = deal(cell(numel(hemiToLoad),1));
for iHemi = 1:numel(hemiToLoad)
    % Loading data for LME analysis
    [gimg,aloc,conds,roiNames] = lme_getData('hemisphere',hemiToLoad{iHemi});
    nconds = numel(conds);
    nlocs = numel(aloc);
    nsubjects = size(gimg,3);
    
    % Generating model predictions
    nNeurons = 360;
    shift = 2.3; % Derived from behavioural data
    % Decision model parameters
    spanDec = 2;
    sigmaDec = 10;
    muDec = rand(nNeurons,1)*spanDec-(spanDec/2);
    funDec = @(m,s) abs(abs(normcdf(aloc,m+s,sigmaDec)-0.5)-0.5)*2;
    
    % Hemifield model parameters
    sigmaHemi = 64; % Based on Salminen et al. 2009
    popCenter = 90;
    propContra = 0.7;
    propIpsi = 1-propContra;
    span = 20;
    % Mean of tuning curves uniformly distributed around +90 and -90 degrees
    % Modelling left hemisphere response here
    if ismember(hemiToLoad{iHemi},{'lh_and_rh','lh'})
        % Contralateral hemifield is Right
        muHemi_contra = rand(round(nNeurons*propContra),1)*span-(span/2)+popCenter;
        muHemi_ipsi = rand(round(nNeurons*propIpsi),1)*span-(span/2)-popCenter;
    else
        % Contralateral hemifield is left
        muHemi_contra = rand(round(nNeurons*propContra),1)*span-(span/2)-popCenter;
        muHemi_ipsi = rand(round(nNeurons*propIpsi),1)*span-(span/2)+popCenter;
    end
    funHemi = @(m,s) normpdf(aloc,m+s,sigmaHemi)/normpdf(m,m,sigmaHemi);
    % Regressors
    % Decisional uncertainty (center-periphery) model
    y = arrayfun(funDec,muDec,zeros(size(muDec)),'UniformOutput',false);
    x_Dec_pre = mean(cat(1,y{:}))';
    y = arrayfun(funDec,muDec,-shift*ones(size(muDec)),'UniformOutput',false);
    x_Dec_postAV = mean(cat(1,y{:}))';
    y = arrayfun(funDec,muDec,shift*ones(size(muDec)),'UniformOutput',false);
    x_Dec_postVA = mean(cat(1,y{:}))';
    
    % Hemifield (left-right) model
    y = cat(1,arrayfun(funHemi,muHemi_contra,zeros(size(muHemi_contra)),'UniformOutput',false),...
        arrayfun(funHemi,muHemi_ipsi,zeros(size(muHemi_ipsi)),'UniformOutput',false));
    x_Sp_pre = mean(cat(1,y{:}))'; % pre-recal left-right contrast
    y = cat(1,arrayfun(funHemi,muHemi_contra,-shift*ones(size(muHemi_contra)),'UniformOutput',false),...
        arrayfun(funHemi,muHemi_ipsi,-shift*ones(size(muHemi_ipsi)),'UniformOutput',false));
    x_Sp_postAV = mean(cat(1,y{:}))';
    y = cat(1,arrayfun(funHemi,muHemi_contra,shift*ones(size(muHemi_contra)),'UniformOutput',false),...
        arrayfun(funHemi,muHemi_ipsi,shift*ones(size(muHemi_ipsi)),'UniformOutput',false));
    x_Sp_postVA = mean(cat(1,y{:}))';
    
    % Psychometric funtion model features
    % Computing group average PMF values
    fLoaded = load('fMRI_behav_all.mat');
    pre = cat(2, fLoaded.out.pre);
    post = cat(2, fLoaded.out.post);
    y = arrayfun(@(i) mean(pre(i).NumPos./pre(i).OutOfNum), 1:size(pre,2),...
        'UniformOutput', false);
    x_Pmf_pre = mean(cat(1,y{:}))';
    % post-AV adaptation
    y = arrayfun(@(i) mean(post(i).NumPos(1:2,:)./post(i).OutOfNum(1:2,:)),...
        1:size(post,2),'UniformOutput', false);
    x_Pmf_postAV = mean(cat(1,y{:}))';
    % post-VA adaptation
    y = arrayfun(@(i) mean(post(i).NumPos(3:4,:)./post(i).OutOfNum(3:4,:)),...
        1:size(post,2),'UniformOutput', false);
    x_Pmf_postVA = mean(cat(1,y{:}))';
    
    nROIs = 1;
    x_Dec_recal = repmat([x_Dec_pre; x_Dec_postAV; x_Dec_postVA], nsubjects*nROIs, 1);
    x_Dec = repmat([x_Dec_pre; x_Dec_pre; x_Dec_pre], nsubjects*nROIs, 1);
    x_Sp_recal = repmat([x_Sp_pre; x_Sp_postAV; x_Sp_postVA], nsubjects*nROIs, 1);
    x_Sp = repmat([x_Sp_pre; x_Sp_pre; x_Sp_pre], nsubjects*nROIs, 1);
    x_Pmf_recal = repmat([x_Pmf_pre; x_Pmf_postAV; x_Pmf_postVA], nsubjects*nROIs, 1);
    x_Pmf = repmat([x_Pmf_pre; x_Pmf_pre; x_Pmf_pre], nsubjects*nROIs, 1);
    
    % Generating lookup table for predictiors
    name = {'Const';'SpRecal';'Sp';'DecRecal';'Dec';'PmfRecal';'Pmf'};
    X = {ones(nsubjects*nlocs*nconds*nROIs,1);x_Sp_recal;x_Sp;...
        x_Dec_recal;x_Dec;x_Pmf_recal;x_Pmf};
    predictorLUT = table(name,X);
    
    model = [];
    % Check if the frst model is 'Const'
    assert(strcmp(modelsToUse{1},'Const'),'The frst model must be ''Const''');
    if numel(modelsToUse) == 3
        modelNames = {modelsToUse{1},...
            sprintf('%s+%s',modelsToUse{2},modelsToUse{1}),...
            sprintf('%s+%s',modelsToUse{3},modelsToUse{1}),...
            sprintf('%s+%s+%s',modelsToUse{2},modelsToUse{3},modelsToUse{1}),...
            sprintf('%sRecal+%s+%s',modelsToUse{2},modelsToUse{3},modelsToUse{1}),...
            sprintf('%s+%sRecal+%s',modelsToUse{2},modelsToUse{3},modelsToUse{1}),...
            sprintf('%sRecal+%sRecal+%s',modelsToUse{2},modelsToUse{3},modelsToUse{1})};
        plotOrder = 4:7;
    elseif numel(modelsToUse) == 4
        modelNames = {modelsToUse{1},...
            sprintf('%s+%s+%s+%s',modelsToUse{2},modelsToUse{3},modelsToUse{4},modelsToUse{1}),...
            sprintf('%sRecal+%s+%s+%s',modelsToUse{2},modelsToUse{3},modelsToUse{4},modelsToUse{1}),...
            sprintf('%s+%sRecal+%s+%s',modelsToUse{2},modelsToUse{3},modelsToUse{4},modelsToUse{1}),...
            sprintf('%s+%s+%sRecal+%s',modelsToUse{2},modelsToUse{3},modelsToUse{4},modelsToUse{1}),...
            sprintf('%s+%sRecal+%sRecal+%s',modelsToUse{2},modelsToUse{3},modelsToUse{4},modelsToUse{1}),...
            sprintf('%sRecal+%s+%sRecal+%s',modelsToUse{2},modelsToUse{3},modelsToUse{4},modelsToUse{1}),...
            sprintf('%sRecal+%sRecal+%s+%s',modelsToUse{2},modelsToUse{3},modelsToUse{4},modelsToUse{1}),...
            sprintf('%sRecal+%sRecal+%sRecal+%s',modelsToUse{2},modelsToUse{3},modelsToUse{4},modelsToUse{1})};
        plotOrder = 2:9;
    end
    nModel = numel(modelNames);
    for iROI=1:5
        for m = 1:nModel
            y = reshape(gimg(:,:,:,iROI), [], 1);
            G = {nominal(reshape(repmat(1:nsubjects, nlocs*nconds, 1, nROIs), [], 1))};
            Z = {ones(nsubjects*nlocs*nconds*nROIs, 1)};
            actModelName = modelNames{m};
            actModelPredictors = strsplit(actModelName,'+');
            % Using a for loop here ensures that the order of models will be
            % the same as in actModelPredictors.
            temp = {};
            for iPred = 1:numel(actModelPredictors)
                temp(iPred) = predictorLUT.X(ismember(predictorLUT.name,...
                    actModelPredictors{iPred}));
            end
            model(m,iROI).X = cat(2,temp{:});
            model(m,iROI).num_par = size(model(m,iROI).X, 2);
            if plotModelDesign
                if iROI == 1
                    figure; imagesc(model(m,iROI).X); colorbar;
                end
            end
            lme = fitlmematrix(model(m,iROI).X, y, Z, G, 'FitMethod', 'ML', ...
                'FixedEffectPredictors', actModelPredictors,...
                'RandomEffectPredictors', {{'Intercept'}}, ...
                'RandomEffectGroups', {'Subject'});
            
            model(m,iROI).lme = lme;
            model(m,iROI).BIC = lme.ModelCriterion.BIC;  % the better the model, the smaller BIC
            model(m,iROI).AIC = lme.ModelCriterion.AIC; % the better the model, the smaller AIC
            model(m,iROI).LL = lme.ModelCriterion.LogLikelihood; % loglike
            model(m,iROI).par = lme.Coefficients.Estimate; % coefficient estimates
            model(m,iROI).pval = lme.Coefficients.pValue; % coefficient p value
        end
    end
    temp = cat(1,model.BIC);
    BIC{iHemi} = cell2table(num2cell(reshape(temp,nModel,5)),...
                    'VariableNames',roiNames,...
                    'RowNames',modelNames(1:nModel));
    temp = cat(1,model.AIC);
    AIC{iHemi} = cell2table(num2cell(reshape(temp,nModel,5)),...
                    'VariableNames',roiNames,...
                    'RowNames',modelNames(1:nModel));
    temp = cat(1,model.LL);
    LL{iHemi} = cell2table(num2cell(reshape(temp,nModel,5)),...
                    'VariableNames',roiNames,...
                    'RowNames',modelNames(1:nModel));
end

if numel(BIC) > 1
    % Averaging BIC, AIC and LL across hemispheres if necessary
    temp = cellfun(@(c) c{:,:},BIC,'UniformOutput',false);
    temp = mean(cat(3,temp{:}),3);
    BIC = BIC{1};
    BIC{:,:} = temp;
    temp = cellfun(@(c) c{:,:},AIC,'UniformOutput',false);
    temp = mean(cat(3,temp{:}),3);
    AIC = AIC{1};
    AIC{:,:} = temp;
    temp = cellfun(@(c) c{:,:},LL,'UniformOutput',false);
    temp = mean(cat(3,temp{:}),3);
    LL = LL{1};
    LL{:,:} = temp;
else
    BIC = BIC{:};
    AIC = AIC{:};
    LL = LL{:};
end
% Here I convert the BIC values computed by Matlab to BIC values according
% to the Bishop definition (this way the relative BIC will be higher if
% the model performs better) 
% see https://uk.mathworks.com/help/econ/aicbic.html for the MATLAB
% definition and Bishop's Pattern Recognition and Machine Learning (2006) 
% p217, eq. 4.139. The scaling factor is -2 between the two. 
temp = BIC{:,:}/-2;
temp = temp-repmat(temp(1,:),size(temp,1),1);
BIC_rel = BIC;
BIC_rel{:,:} = temp;

out = struct;
out.BIC = BIC;
out.BIC_rel = BIC_rel;
out.AIC = AIC;
out.LL = LL;

if plotFigure
    % Writing actual BIC values on top of bars
    pltYvalTxt = true;
    % plot relative BIC
    nROIs = numel(roiNames);
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
    figure();
    set(gcf, 'Position', [0 0 1000 200], 'PaperPositionMode', 'auto');
    temp = BIC_rel{:,:};
    yl = [min(temp(:)),max(temp(:))]*1.05;
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
    temp_model = repmat(modelNames(plotOrder)',numel(roiNames),1);
    temp_data = BIC_rel{plotOrder,:};
    logBF_table = table(temp_roi,temp_model,temp_data(:),'VariableNames',...
        {'roi','model','logBF'});
    varargout = {out,logBF_table};
else
    varargout = {out};
end

end
