function pf_table = figure_PF_prepost(dataPath)

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

% Fine grained stimulus location
StimLevelsFineGrain=min(StimLevels):max(StimLevels)/1000:max(StimLevels);

% Initialize table
pf_table = table(reshape(repmat({'pre' 'postAV' 'postVA'}, ...
    numel(StimLevelsFineGrain), 1), [], 1), ...
    repmat(StimLevelsFineGrain', 3, 1), ...
    'VariableNames',{'condition','auditory_location'});

% Number of subjects
nsubjects = numel(out);

% Open figure and define colours
fig = figure('Name', modality);
color = {[133,192,249]./255,[169,90,161]./255,[15,32,128]./255};

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
    
    % Initial fit (needed for multifit comparison)
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
    
    % Plot individual PF
    for i=1:3
        plot_PF(fig, PF, mdl_fuller(ss).params(i,:), StimLevels, color{i}, 0);
    end    
    
    % Add individual PF to table
    pf_table.(sprintf('sub%d_percent_perceived_right', ss)) = [...
        PF(mdl_fuller(ss).params(1,:), StimLevelsFineGrain) ...
        PF(mdl_fuller(ss).params(2,:), StimLevelsFineGrain) ...
        PF(mdl_fuller(ss).params(3,:), StimLevelsFineGrain)]';
end

% Calculate group level parameters
params = cat(3, mdl_fuller.params);

% Plot group PF
for i=1:3
    plot_PF(fig, PF, mean(params(i,:,:), 3), StimLevels, color{i}, 1);
end
legend({'pre' 'postAV' 'postVA'}, 'Location', 'NorthWest');
title(modality);

% Add group PF to table
pf_table.mean_percent_perceived_right = [...
    PF(mean(params(1,:,:), 3), StimLevelsFineGrain) ...
    PF(mean(params(2,:,:), 3), StimLevelsFineGrain) ...
    PF(mean(params(3,:,:), 3), StimLevelsFineGrain)]';


function plot_PF(fig, PF, paramsValues, StimLevels, color, flag)

figure(fig);
set(fig, 'PaperPositionMode', 'auto');
hold on;

% Calculate PF at fine resolution
StimLevelsFineGrain=min(StimLevels):max(StimLevels)/1000:max(StimLevels);
ProportionCorrectModel = PF(paramsValues(1,:), StimLevelsFineGrain);

% Plot PF
if flag
    plot(StimLevelsFineGrain, ProportionCorrectModel, 'Color', color, ...
        'LineWidth', 4);
else
    plot(StimLevelsFineGrain, ProportionCorrectModel, 'Color', color, ...
        'LineWidth', 0.5);
end

% Update axes, legend etc
set(gca, 'FontSize', 16, 'Xtick', StimLevels, 'YTick', [0:0.2:1], 'box', ...
    'off');
xlabel(gca, 'A location (°)', 'FontSize', 20);
ylabel(gca, '% A decoded right', 'FontSize', 20);
xlim(gca, [-15 15]);
ylim(gca, [0 1]);
pbaspect([1 1 1]);