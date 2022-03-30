function [data, scaleparams] = mvpa_norm_calc(data, cfg, scaleparams)
%   mvpa_norm_calc(data, cfg, scaleparams)
% see also decoding_scale_data in TDT

if ~exist('scaleparams', 'var')
    switch cfg.scale.method
        case 'min0max1'
            scaleparams.samples_min = min(data, [], 1);
            scaleparams.samples_max = max(data, [], 1);
            min_eq_max = scaleparams.samples_min==scaleparams.samples_max; % check if in any dimension min == max [taken from TDT]
            scaleparams.samples_max(min_eq_max) = scaleparams.samples_min(min_eq_max) + 1; % prevents divide by 0, if min == max [taken from TDT]
        case 'mean'
            scaleparams.samples_mean = mean(data);
        case 'z'
            scaleparams.samples_mean = mean(data);
            scaleparams.samples_std = std(data);
            scaleparams.samples_std(scaleparams.samples_std==0) = 1; % prevents divide by 0, if no std exists [taken from TDT]
    end
end

% Feature mean centering or z-transformation
if ismember(cfg.scale.method, {'z' 'mean'})
    for i=1:size(data, 1)
        data(i,:) = data(i,:) - scaleparams.samples_mean;
        if strcmp(cfg.scale.method, 'z')
            data(i,:) = data(i,:) ./ scaleparams.samples_std;
        end
    end
    data(data<cfg.scale.cutoff(1)) = cfg.scale.cutoff(1);
    data(data>cfg.scale.cutoff(2)) = cfg.scale.cutoff(2);
end

% Feature scaling into [0 1] range
if strcmp(cfg.scale.estimation, 'min0max1')
    for i=1:size(data, 1)
        data(i,:) = (data(i,:) - scaleparams.samples_min) ./ (scaleparams.samples_max - scaleparams.samples_min);
    end
end