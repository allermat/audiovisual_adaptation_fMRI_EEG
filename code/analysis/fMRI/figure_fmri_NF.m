function [fmri_NF_table, fmri_NF_data_table] = figure_fmri_NF(dataPath)

% Psychometric function
PF = @PAL_CumulativeNormal;

% Load input data
load(dataPath);

nsubjects = numel(out.data);
nROIs = numel(out.ROI);

% Stimulus locations
StimLevels = [-12 -5 -2 0 2 5 12];

% Fine grained stimulus location
StimLevelsFineGrain=min(StimLevels):max(StimLevels)/1000:max(StimLevels);

% Initialize tables
fmri_NF_table = table(reshape(repmat(out.ROI', ...
    numel(StimLevelsFineGrain) * 3, 1), [], 1), ...
    reshape(repmat({'pre' 'postAV' 'postVA'}, numel(StimLevelsFineGrain), ...
    nROIs), [], 1), reshape(repmat(StimLevelsFineGrain', 3, nROIs), [], 1), ...
    'VariableNames',{'roi','condition','auditory_location'});
fmri_NF_data_table = table(reshape(repmat(out.ROI', ...
    numel(StimLevels) * 3, 1), [], 1), ...
    reshape(repmat({'pre' 'postAV' 'postVA'}, numel(StimLevels), ...
    nROIs), [], 1), reshape(repmat(StimLevels', 3, nROIs), [], 1), ...
    'VariableNames',{'roi','condition','auditory_location'});

% Format data
for i=1:nROIs
    for ss=1:nsubjects
        for j=1:3
            OutOfNum{i}(j,:,ss) = out.data(ss).OutOfNum{i,j};
            NumPos{i}(j,:,ss) = out.data(ss).NumPos{i,j};
        end
    end
    OutOfNum{i} = sum(OutOfNum{i}, 3);
    NumPos{i} = sum(NumPos{i}, 3);
end

color = {[133,192,249]./255,[169,90,161]./255,[15,32,128]./255};

for i=1:nROIs    
    % Initial fit to get parameter estimates
    for j=1:3
        paramsInit(j,:) = fit_PF(PF, StimLevels, NumPos{i}(j,:), ...
            OutOfNum{i}(j,:), [], get_PAL_opt_precise);
    end
    
    % Eta parameter for betabinomial model
    eta = 0.15;

    % Fit fuller model (recalibration model)
    thetasInit = [paramsInit(:,1)' mean(paramsInit(:,2)) paramsInit(1,4)];
    paramsIDmatrix = [1 4 5 5; 2 4 5 5; 3 4 5 5];
    mdl_fuller(i) = multifit_PF(PF, 'fuller model', thetasInit, ...
        paramsIDmatrix, StimLevels, NumPos{i}, OutOfNum{i}, [0 0.45], eta);
    
    fig = figure('Name', sprintf('group %s', out.ROI{i}));
    for j=1:3
        plot_PF(fig, PF, mdl_fuller(i).params(j,:), StimLevels, ...
            out.ROI{i}, color{j});
        plot_data(fig, StimLevels, NumPos{i}(j,:), OutOfNum{i}(j,:), color{j})
    end
    
    % Save individual NF to variable
    decoded_right{i} = [PF(mdl_fuller(i).params(1,:), StimLevelsFineGrain) ...
        PF(mdl_fuller(i).params(2,:), StimLevelsFineGrain) ...
        PF(mdl_fuller(i).params(3,:), StimLevelsFineGrain)]';
end

% Add NF to table
fmri_NF_table.percent_decoded_right = cat(1, decoded_right{:});

% Add data points underlying NF fit to table
clear decoded_right
for i=1:nROIs
    for j=1:3
        decoded_right{j,i}(:,1) = NumPos{i}(j,:) ./ OutOfNum{i}(j,:);
    end
end
fmri_NF_data_table.percent_decoded_right = cat(1, decoded_right{:});


function plot_PF(fig, PF, paramsValues, StimLevels, ROI, color)

figure(fig);
set(fig, 'PaperPositionMode', 'auto');

StimLevelsFineGrain=min(StimLevels):max(StimLevels)/1000:max(StimLevels);
ProportionCorrectModel = PF(paramsValues(1,:), StimLevelsFineGrain);
hold on;
plot(StimLevelsFineGrain, ProportionCorrectModel, 'Color', color, 'LineWidth', 2);
set(gca, 'FontSize', 18, 'Xtick', StimLevels, 'YTick', [0:0.2:1], 'box', 'off');
xlabel(gca, 'A location (°)', 'FontSize', 18);
ylabel(gca, '% A decoded right', 'FontSize', 24);
title(ROI, 'FontSize', 18);
xlim(gca, [-15 15]);
ylim(gca, [0 1]);
pbaspect([1 1 1]);


function plot_data(fig, StimLevels, NumPos, OutOfNum, color)

figure(fig);

ProportionCorrectObserved = NumPos ./ OutOfNum;
hold on;
plot(StimLevels, ProportionCorrectObserved, 'ko','MarkerSize', 10, ...
    'MarkerFaceColor', color);
