function pse_table = figure_PSE_prepost(dataPath)

% Psychometric function
PF = @PAL_CumulativeNormal;

% Load input data
load(dataPath);

% Data modality
[pathstr, name, ext] = fileparts(dataPath); 
tmp = strsplit(name, '_');
modality = tmp{1};

% Stimulus locations
StimLevels = [-12 -5 -2 0 2 5 12];

% Initialize table
pse_table = table({'pre'; 'postAV'; 'postVA'}, 'VariableNames',{'condition'});

% Number of subjects
nsubjects = numel(out);

for ss=1:nsubjects
    % Pretest
    NumPos{ss}(1,:) = sum(out(ss).pre.NumPos, 1);
    OutOfNum{ss}(1,:) = sum(out(ss).pre.OutOfNum, 1);
    
    % Posttest
    if ss == 1 && strcmp(modality, 'EEG')
        % postAV adaptation
        NumPos{ss}(2,:) = out(ss).post.NumPos(1,:);
        OutOfNum{ss}(2,:) = out(ss).post.OutOfNum(1,:);
        % postVA adaptation
        NumPos{ss}(3,:) = out(ss).post.NumPos(2,:);
        OutOfNum{ss}(3,:) = out(ss).post.OutOfNum(2,:);
    else
        % postAV adaptation
        NumPos{ss}(2,:) = sum(out(ss).post.NumPos(1:2,:), 1);
        OutOfNum{ss}(2,:) = sum(out(ss).post.OutOfNum(1:2,:), 1);
        % postVA adaptation
        NumPos{ss}(3,:) = sum(out(ss).post.NumPos(3:4,:), 1);
        OutOfNum{ss}(3,:) = sum(out(ss).post.OutOfNum(3:4,:), 1);
    end
    
    % Initial fit
    for i=1:3
        paramsInit(i,:) = fit_PF(PF, StimLevels, NumPos{ss}(i,:), ...
            OutOfNum{ss}(i,:), [], get_PAL_opt_precise);
    end
    
    % Eta parameter for betabinomial model
    eta = 0.15;

    % Fit fuller model (recalibration model)
    thetasInit = [paramsInit(:,1)' mean(paramsInit(:,2)) paramsInit(1,4)];
    paramsIDmatrix = [1 4 5 5; 2 4 5 5; 3 4 5 5];
    mdl_fuller(ss) = multifit_PF(PF, 'fuller model', thetasInit, ...
        paramsIDmatrix, StimLevels, NumPos{ss}, OutOfNum{ss}, [0 0.1], eta);
end

% Calculate group level parameters
params = cat(3, mdl_fuller.params);
pse = squeeze(params(:,1,:));

% Open figure and define colours
fig = figure('Name', modality);
color = {[133,192,249]./255,[169,90,161]./255,[15,32,128]./255};

% Plot PSE estimates
plot_pse(fig, pse, color)
title(modality);

% Add PSE values to table
pse_table.mean_pse = mean(pse, 2);
pse_table.sem_pse = std(pse, [], 2) / sqrt(nsubjects);
for ss=1:nsubjects
    pse_table.(sprintf('sub%d_pse', ss)) = pse(:,ss);
end


function plot_pse(fig, data, color)

figure(fig);
set(fig, 'PaperPositionMode', 'auto');
hold on;

% Calculate mean, SEM
nsubjects = size(data, 2);
mean_data = mean(data, 2);
sem_data = std(data, [], 2) / sqrt(nsubjects);

% Plot bars
oo = [2 1 3];
for o=1:numel(oo)
    bar(oo(o), mean_data(o), 'FaceColor', color{o});
    errorbar(oo(o), mean_data(o), sem_data(o), 'Color', 'k');
    plot(ones(1, nsubjects)*0.8+oo(o)-1, data(o,:), 'ko', 'MarkerSize', 10);
end

% Plot connecting lines
for i=1:nsubjects
    plot(0.8:1:2.8, data(oo,i), 'k--');
end

% Update axes, legend etc
xlim([0 4]);
ylim([-5 5]);
ylabel('PSE', 'FontSize', 20);
set(gca, 'XTick', 1:3, 'XTickLabel', {'postAV' 'pre' 'postVA'}, ...
    'YTick', [-4:2:4], 'FontSize', 16, 'box', 'off');
pbaspect([1 1 1]);