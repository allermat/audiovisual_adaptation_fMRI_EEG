function [paramsValues, lineseries, SD, paramsSimFilt] = fit_PF(PF, StimLevels, NumPos, OutOfNum, opt, PAL_opt)

lineseries = [];
SD = [];

if ~exist('opt', 'var') || isempty(opt)
    opt = struct();
end

% Add palamedes toolbox to the path
if isempty(regexp(path,'Palamedes','once'))
    paldir = fullfile(get_path('toolbox'), 'Palamedes');
    if isdir(paldir)
        addpath(paldir)
    else
        error('Palamedes toolbox is not in the path.')
    end
end

% Parameter grid defining parameter space through which to perform a brute-force
% search for values to be used as initial guesses in iterative parameter search
params = {'alpha' 5:1:35; 'beta' 0.1:0.1:1; 'gamma' 0.02; 'lambda' 0.02};
for i=1:size(params, 1)
    if isfield(opt, params{i,1})
        searchGrid.(params{i,1}) = opt.(params{i,1}); % value defined by opt
    else
        searchGrid.(params{i,1}) = params{i,2}; % default value
    end
end
paramsFree = ~structfun(@isscalar, searchGrid);  %1: free parameter, 0: fixed parameter [1 1 1 1]

% Optional PAL settings
if ~exist('PAL_opt', 'var')
    PAL_opt = PAL_minimize('options');   %type PAL_minimize('options','help') for help
    PAL_opt.TolFun = 1e-09;     %increase required precision on LL
    PAL_opt.MaxIter = 100;
    PAL_opt.Display = 'off';    %suppress fminsearch messages
end

% Perform fit
[paramsValues, LL, exitflag, output] = PAL_PFML_Fit(StimLevels, NumPos, ...
    OutOfNum, searchGrid, paramsFree, PF, 'searchOptions', PAL_opt);