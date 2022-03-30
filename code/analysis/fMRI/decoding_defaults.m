% function cfg = decoding_defaults(cfg)
%
% Function where all defaults are declared and paths are set for decoding
% toolbox.
%
% Usage:
%  cfg = decoding_defaults
%       Get all decoding defaults (OVERWRITES cfg)
%
%  cfg = decoding_defaults(cfg)
%       Puts default values to all fields that are not defined in the
%       passed cfg (UPDATES cfg).
%
% See also DECODING

% Martin H. 2011/03

% HISTORY
% MARTIN, 11/12/18
%   added self-referencing function to add values from defaults to cfg that
%   have not been set manually
% KAI, 11/07/01
%   wrap-around correction and cfg.searchlight.wrap_control added
% A MIHALIK, 16/12/31
%   path and some default values changed
% A MIHALIK, 17/03/16
%   path, ROI mod, image scaling added
% A MIHALIK, 19/03/22
%   TDT not added to path anymore

function cfg = decoding_defaults(cfg)

if ~exist('cfg','var'), cfg = struct; end

%% Add paths
defaults.report = []; % init field

% path to libSVM (set if version delivered with this package is not working)
% addpath('/analysis/share/software/matlab_libraries/libsvm-3.11')

%% Set defaults

% General values
defaults.testmode = 0; % Test mode off
% defaults.analysis = 'searchlight'; % 'searchlight' 'ROI', 'wholebrain'
defaults.software = 'SPM12';

% display options
defaults.verbose = -1; % Verbosity (0 to 2)
defaults.plot_design = 1; % decide whether you want to save the design as image. 
                                %     0: no plotting (not recommended)
                                %     1: plot using the default files formats
                                %     2: will be plotted only at the end
% default.plot_design_formats = {'-dpng', '-depsc2'}; % list all formats that you want to save the figure as
defaults.plot_selected_voxels = 1; % a value of n means that the currently 
                                % selected voxels (e.g. a searchlight, ROI, 
                                % ...) are plotted every n-th step 

% specification of scaling (TDT extended with image scaling and feature centering!)
defaults.scale.method = 'none'; % none z min0max1 mean [PRONTO default]
defaults.scale.estimation = 'none'; % none all across separate
defaults.scale.cutoff = [-inf inf]; % -inf inf -3 3
defaults.scale.image = 'eucledian'; % none z eucledian mean [PRONTO default]

% Specification of feature transformation
defaults.feature_transformation.method = 'none';
defaults.feature_transformation.estimation = 'none';

% Specification of feature selection
defaults.feature_selection.method = 'none';
defaults.feature_selection.estimation = 'across'; % i.e. carry out only on training data and apply to test data
defaults.feature_selection.optimization_criterion = 'max';

% Specification of parameter selection
% defaults.parameter_selection.design.function.name = 'make_design_cv_misc';
defaults.parameter_selection.method = 'none'; % none grid
% defaults.parameter_selection.parameters = {'-c' '-n'};
% defaults.parameter_selection.parameter_range = {[0.0001 0.001 0.01 0.1 1 10 100 10000]; [0.001 0.01 0.1:0.1:1]};
defaults.parameter_selection.format.name = 'string_number';
defaults.parameter_selection.format.separator = ' ';
defaults.parameter_selection.optimization_criterion = 'min';

% Searchlight specific defaults
defaults.searchlight.unit = 'mm'; % searchlight unit ('mm' or 'voxels')
defaults.searchlight.radius = 12; % 4 voxels is a standard often used
defaults.searchlight.spherical = 0;
defaults.searchlight.wrap_control = 1; % tests that no wrap-around effects occur when searchlight is shifted.
% Only switch off when you are sure that your brain is not near the border and when you need even more speed (the check is extremely fast, though).

% ROI specific defaults
defaults.roi.extension = '.nii';
defaults.roi.mask.atlas = repmat([repmat({'FS_Destrieux'}, 3, 1); ...
    repmat({fullfile('ProbAtlas_v4', 'freesurfer10')}, 2, 1)], 3, 1);
defaults.roi.mask.name = {'HG'; 'hA'; 'IPL'; 'IPS'; 'FEF'}; 
defaults.roi.mask.name = [defaults.roi.mask.name; ...
    cellfun(@(x) ['lh-' x], defaults.roi.mask.name, 'un', 0); ...
    cellfun(@(x) ['rh-' x], defaults.roi.mask.name, 'un', 0)];
defaults.roi.mask.fname = {'G_temp_sup-G_T_transv'; 'hAud'; 'IPL'; ...
    'IPS0-5_maxprob'; 'FEF_maxprob'};
defaults.roi.mask.fname = [cellfun(@(x) ['rSROI_' x], ...
    defaults.roi.mask.fname, 'un', 0); cellfun(@(x) ['rSROI_lh.' x], ...
    defaults.roi.mask.fname, 'un', 0); cellfun(@(x) ['rSROI_rh.' x], ...
    defaults.roi.mask.fname, 'un', 0)];
if isfield(cfg, 'roi') && ~isfield(cfg.roi, 'fileindex') || ~isfield(cfg, 'roi')
    cfg.roi.fileindex = 1:numel(defaults.roi.mask.fname);
end
tmp = {'atlas' 'name' 'fname'};
for i=1:numel(tmp)
    defaults.roi.mask.(tmp{i}) = defaults.roi.mask.(tmp{i})(cfg.roi.fileindex);
end

% Decoding specific values for libsvm
defaults.decoding.software = 'libsvm'; % libsvm as a standard
defaults.decoding.kernel.function = @(X,Y) X*Y'; % for kernel method linear kernel as default
defaults.decoding.kernel.pass_vectors = 0; % if 1, original data vectors will be passed
                                           % in addition to the kernel as data_train./_test.vectors .
                                           % might be useful if you e.g. need the dimension of the 
                                           % original data
if isfield(cfg, 'decoding') && isfield(cfg.decoding, 'machine')
   switch cfg.decoding.machine
       case 'C-SVC'
           defaults.decoding.train.classification.model_parameters = '-s 0 -t 0 -c 1 -b 0 -q'; % linear classification with C-SVC
       case 'nu-SVC'
           defaults.decoding.train.classification.model_parameters = '-s 0 -t 0 -n 0.5 -b 0 -q'; % linear classification with nu-SVC
       case 'nu-SVR'
            defaults.decoding.train.regression.model_parameters = '-s 4 -t 0 -c 1 -n 0.5 -b 0 -q'; % nu-SVR (adapt cost to control speed)
       case 'e-SVR'
           defaults.decoding.train.regression.model_parameters = '-s 3 -t 0 -c 1 -p 0.1 -b 0 -q'; % nu-SVR (adapt cost to control speed)
   end
else
    defaults.decoding.train.classification_kernel.model_parameters = '-s 0 -t 4 -c 1 -b 0 -q'; % linear classification with kernel
end
defaults.decoding.method = char(fieldnames(defaults.decoding.train));
defaults.decoding.test.classification.model_parameters = '-q';
defaults.decoding.test.classification_kernel.model_parameters = '-q'; % linear classification
defaults.decoding.test.regression.model_parameters = '-q';

% Results specific defaults
defaults.results.output = {'predicted_labels'};
defaults.results.write = 1; % 1: write results both as .mat and as image, 2: mat only, 0: dont write
defaults.results.backgroundvalue = 0; % background of images consists of zeros
defaults.results.overwrite = 0; % overwrite existing results!!
defaults.results.setwise = 1; % return results of each decoding set separately

%% Add values to cfg that have not yet been set
cfg = assign_fields(defaults,cfg);

%==============================================
function cfg = assign_fields(defaults,cfg)

% Self-referencing function that goes through all field names and adds 
% non-existent fields to cfg from the defaults.

d_fields = fieldnames(defaults);

for i = 1:size(d_fields,1)
    % If there are no subfields in the current field
    if ~isstruct(defaults.(d_fields{i}))
        % If this field doesn't exist in cfg, add it from defaults
        if ~isfield(cfg,d_fields{i})
            cfg.(d_fields{i}) = defaults.(d_fields{i});
        end
    % If there are subfields in the current field
    else
        % If this field doesn't exist in cfg, add it (and all subfields) from defaults
        if ~isfield(cfg,d_fields{i})
            cfg.(d_fields{i}) = defaults.(d_fields{i});
        % Else loop through function again for all subfields    
        else
            cfg.(d_fields{i}) = assign_fields(defaults.(d_fields{i}),cfg.(d_fields{i}));
        end
    end
end