function varargout = run_lme_fMRI_pre(varargin)

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
    % Only using pre-test
    nconds = 1;
    nlocs = numel(aloc);
    nsubjects = size(gimg,3);
    % Model parameters
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
    
    % Hemifield (left-right) model
    y = cat(1,arrayfun(funHemi,muHemi_contra,zeros(size(muHemi_contra)),'UniformOutput',false),...
        arrayfun(funHemi,muHemi_ipsi,zeros(size(muHemi_ipsi)),'UniformOutput',false));
    x_Sp_pre = mean(cat(1,y{:}))'; % pre-recal left-right contrast
    
    % Psychometric funtion model features
    % Computing group average PMF values
    fLoaded = load('fMRI_behav_all.mat');
    pre = cat(2, fLoaded.out.pre);
    y = arrayfun(@(i) mean(pre(i).NumPos./pre(i).OutOfNum), 1:size(pre,2),...
        'UniformOutput', false);
    x_Pmf_pre = mean(cat(1,y{:}))';
    
    nROIs = 1;
    x_Dec_noshift = repmat(x_Dec_pre, nsubjects*nROIs, 1);
    x_Sp_noshift = repmat(x_Sp_pre, nsubjects*nROIs, 1);
    x_Pmf_noshift = repmat(x_Pmf_pre, nsubjects*nROIs, 1);
    
    % Generating lookup table for predictiors
    name = {'Const';'Sp';'Dec';'Pmf'};
    X = {ones(nsubjects*nlocs*nconds*nROIs,1);x_Sp_noshift;...
        x_Dec_noshift;x_Pmf_noshift};
    predictorLUT = table(name,X);
    
    % Check if the frst model is 'Const'
    assert(strcmp(modelsToUse{1},'Const'),'The frst model must be ''Const''');
    if numel(modelsToUse) == 3
        modelNames = {modelsToUse{1},...
            sprintf('%s+%s',modelsToUse{2},modelsToUse{1}),...
            sprintf('%s+%s',modelsToUse{3},modelsToUse{1}),...
            sprintf('%s+%s+%s',modelsToUse{2},modelsToUse{3},modelsToUse{1})};
        plotOrder = {2:4};
    elseif numel(modelsToUse) == 4
        modelNames = {modelsToUse{1},...
            sprintf('%s+%s',modelsToUse{2},modelsToUse{1}),...
            sprintf('%s+%s',modelsToUse{3},modelsToUse{1}),...
            sprintf('%s+%s+%s',modelsToUse{2},modelsToUse{3},modelsToUse{1}),...
            sprintf('%s+%s',modelsToUse{4},modelsToUse{1}),...
            sprintf('%s+%s+%s',modelsToUse{2},modelsToUse{4},modelsToUse{1}),...
            sprintf('%s+%s+%s',modelsToUse{3},modelsToUse{4},modelsToUse{1}),...
            sprintf('%s+%s+%s+%s',modelsToUse{2},modelsToUse{3},modelsToUse{4},modelsToUse{1})};
        plotOrder = {2:8};
    end
    
    nModel = numel(modelNames);
    model = {};
    for iROI=1:numel(roiNames)
        for m = 1:nModel
            y = reshape(gimg(:,1,:,iROI), [], 1);
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
            model{iROI}(m).X = cat(2,temp{:});
            model{iROI}(m).num_par = size(model{iROI}(m).X, 2);
            if plotModelDesign
                if iROI == 1
                    figure; imagesc(model{iROI}(m).X); colorbar;
                end
            end
            lme = fitlmematrix(model{iROI}(m).X, y, Z, G, 'FitMethod', 'ML', ...
                'FixedEffectPredictors', actModelPredictors,...
                'RandomEffectPredictors', {{'Intercept'}}, ...
                'RandomEffectGroups', {'Subject'});
            
            model{iROI}(m).lme = lme;
            model{iROI}(m).BIC = lme.ModelCriterion.BIC;  % the better the model, the smaller BIC
            model{iROI}(m).AIC = lme.ModelCriterion.AIC; % the better the model, the smaller AIC
            model{iROI}(m).LL = lme.ModelCriterion.LogLikelihood; % loglike
            model{iROI}(m).par = lme.Coefficients.Estimate; % coefficient estimates
            model{iROI}(m).pval = lme.Coefficients.pValue; % coefficient p value
        end
    end
    temp = cellfun(@(c) cat(1,c.BIC),model,'UniformOutput',false);
    BIC{iHemi} = array2table(cat(2,temp{:}),'VariableNames',roiNames,...
                    'RowNames',modelNames(1:nModel));
    temp = cellfun(@(c) cat(1,c.AIC),model,'UniformOutput',false);
    AIC{iHemi} = array2table(cat(2,temp{:}),'VariableNames',roiNames,...
                    'RowNames',modelNames(1:nModel));
    temp = cellfun(@(c) cat(1,c.LL),model,'UniformOutput',false);
    LL{iHemi} = array2table(cat(2,temp{:}),'VariableNames',roiNames,...
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
    nROIs = numel(roiNames);
    % plot BIC
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
    % Writing actual BIC values on top of bars
    pltYvalTxt = true;
    figure();
    nRows = numel(plotOrder);
    set(gcf, 'Position', [0 0 1000 200], 'PaperPositionMode', 'auto');
    temp = BIC_rel{:,:};
    yl = [min(temp(:)),max(temp(:))]*1.05;
    for iRow = 1:nRows
        for iROI = 1:nROIs
            subplot(nRows,nROIs,(iRow-1)*nROIs+iROI);
            for k = 1:numel(plotOrder{iRow})
                bar(k,BIC_rel{plotOrder{iRow}(k),iROI},'FaceColor',colors{k}); hold on;
            end
            if pltYvalTxt
                Y = BIC_rel{plotOrder{iRow},iROI};
                ycoord = Y;
                ycoord(ycoord < 0) = 0;
                text(1:length(Y),ycoord,num2str(Y,'%0.1f'),'HorizontalAlignment','center',...
                    'VerticalAlignment','bottom','FontSize',7,'Margin',0.01,...
                    'BackgroundColor','w');
            end
            if iROI == 1, ylabel('logBF'); end
            ylim(yl);
            set(gca,'XTickLabels',[]);
            if iRow == 1, title(roiNames{iROI}); end
        end
        legend(BIC_rel.Properties.RowNames(plotOrder{iRow}));
    end
    temp_roi = repmat(roiNames,numel(modelNames),1);
    temp_roi = temp_roi(:);
    temp_model = repmat(modelNames',numel(roiNames),1);
    temp_data = BIC_rel{:,:};
    logBF_table = table(temp_roi,temp_model,temp_data(:),'VariableNames',...
        {'roi','model','logBF'});
    varargout = {out,logBF_table};
else
    varargout = {out};
end

end
