function [fmri_mvpa_Acc_table, fmri_mvpa_RI_table] = figure_fmri_mvpa_index(dataPath)

% Load data
load(dataPath);

nsubjects = numel(out.data);
nROIs = numel(out.ROI);

% Get pretest and posttest indexes
accuracy = cat(2, out.data.accuracy);
recal_index = squeeze(mean(cat(3, out.data.recal_index), 2));

% Plot data
for i=1:nROIs
    fig = figure('Name', sprintf('group %s', out.ROI{i}));
    plot_estimates(fig, accuracy(i,:), recal_index(i,:), out.ROI{i}, ...
        [-0.4 1.4; -0.2 0.7])    
end

% Initialize output tablew
fmri_mvpa_Acc_table = table(out.ROI, 'VariableNames',{'ROI'});
fmri_mvpa_RI_table = table(out.ROI, 'VariableNames',{'ROI'});

% Add data to accuracy table
fmri_mvpa_Acc_table.mean_accuracy = mean(accuracy, 2);
fmri_mvpa_Acc_table.sem_accuracy = std(accuracy, [], 2) / sqrt(nsubjects);
for j=1:nsubjects
    fmri_mvpa_Acc_table.(sprintf('sub%d_accuracy', j)) = accuracy(:,j);
end

% Add data to recalibration index table
fmri_mvpa_RI_table.mean_recalibration_index = mean(recal_index, 2);
fmri_mvpa_RI_table.sem_recalibration_index = std(recal_index, [], 2) / sqrt(nsubjects);
for j=1:nsubjects
    fmri_mvpa_RI_table.(sprintf('sub%d_recalibration_index', j)) = ...
        recal_index(:,j);
end

function plot_estimates(fig, dataLeft, dataRight, ROI, ylimit)

figure(fig);
set(fig, 'PaperPositionMode', 'auto');

nsubjects = size(dataLeft, 2);

mean_left = mean(dataLeft, 2);
sem_left = std(dataLeft, [], 2) / sqrt(nsubjects);
hold on;
bar(1, mean_left, 'FaceColor', [0.8 0.8 0.8]); % 'none'
errorbar(1, mean_left, sem_left, 'Color', 'k');
plot(ones(1, nsubjects)*0.8, dataLeft, 'ko', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 10);
xlim([0 4]);
ylim(ylimit(1,:));
ylabel('R (Fisher transformed)', 'FontSize', 18);
set(gca, 'FontSize', 18, 'XTick', [], 'box', 'off', 'YTick', [0 0.5 1]);
pbaspect([1 1 1]);

mean_right = mean(dataRight, 2);
sem_right = std(dataRight, [], 2) / sqrt(nsubjects);
pos = get(gca, 'Position');
ax2 = axes('Position', pos, 'XAxisLocation', 'bottom', 'YAxisLocation', ...
    'right', 'box', 'off', 'Color', 'none');
hold on;
bar(3, mean_right, 'FaceColor', [0.8 0.8 0.8]); % 'none'
errorbar(3, mean_right, sem_right, 'Color', 'k');
plot(ones(1, nsubjects)*0.8+2, dataRight, 'ko', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 10);
xlim([0 4]);
ylim(ylimit(2,:));
ylabel('RI', 'FontSize', 18); % Recalibration index
set(ax2, 'FontSize', 18, 'XTick', [1 3], 'XTickLabel', ...
    {'Space' 'Recalibration'}, 'YTick', [0 0.2 0.4]);
pbaspect([1 1 1]);

title(gca, ROI, 'FontSize', 18);