function files = mvpa_make_files(cfg, LUT)
%   mvpa_make_files prepares cfg.files structure with duplicating data to enable
% separate filtering for training and test data

% ------------- Build up initial files structure -------------

% --- 1. cfg.files.name: a nx1 cell array of file names
% files.name = LUT.beta;

% --- 2. cfg.files.chunk: a nx1 vector to indicate what data you want to keep
%       together for cross-validation (typically runs, so enter run numbers)
files.chunk = LUT.chunk;
% files.chunk = arrayfun(@(x) find(unique(LUT.session)==x), LUT.session);

% --- 3. cfg.files.label %: a nx1 vector of labels (for decoding, you can choose 
%       any two numbers as class labels, but normally we use 1 and -1)
if strfind(cfg.decoding.method, 'classification')
    files.label = sign(LUT.aloc);
elseif strcmp(cfg.decoding.method, 'regression') && ~isempty(strfind(cfg.results.dir, 'abs'))
    files.label = abs(LUT.aloc);
else
    files.label = LUT.aloc;
end

% --- 4. cfg.files.set
files.set = ones(size(files.chunk));

% Additional file fields
if strfind(cfg.decoding.method, 'classification')
    files.aloc = LUT.aloc;
end
% if cfg.decoding.scheme == 3
files.blocktype = LUT.blocktype;
% end

% ------------ Filter files by example selection ------------

% Example selection
set_train_filter = {};
set_test_filter = {};
if ismember('aloc5.12', cfg.example_selection)
    set_train_filter{end+1} = ismember(LUT.aloc, [-12 -5 5 12]);
end
if ismember('alocL', cfg.example_selection)
    set_train_filter{end+1} = LUT.aloc < 0;
elseif ismember('alocR', cfg.example_selection)
    set_train_filter{end+1} = LUT.aloc > 0;
elseif ismember('alocM', cfg.example_selection)
    set_train_filter{end+1} = LUT.aloc < 10 & LUT.aloc > -10;
elseif ismember('alocS', cfg.example_selection)
    set_train_filter{end+1}  = LUT.aloc > 10 | LUT.aloc < -10;
end
if strfind(cfg.decoding.method, 'classification')
    [set_train_filter{end+1}, set_test_filter{end+1}] = deal(LUT.aloc ~= 0);
end
set_train_filter{end+1} = ismember(LUT.blocktype, 'pretest');
set_test_filter{end+1} = ismember(LUT.blocktype, {'posttest-radapt' 'posttest-ladapt'}); % & LUT.resp == 0
set_filter = all(cell2mat(set_train_filter), 2) | all(cell2mat(set_test_filter), 2);

% Define example classes for CV and generalization
files.cvclass = zeros(size(files.chunk));
files.cvclass(ismember(LUT.blocktype, 'pretest')) = 1; % all examples in CV loop
files.cvclass(ismember(LUT.blocktype, 'pretest') & LUT.resp == 0) = 2; % examples in CV loop, but with no response
files.cvclass(~ismember(LUT.blocktype, 'pretest') & LUT.resp == 0) = 3; % examples not in CV loop and with no response

% Filter files
% fieldname = setdiff(fieldnames(files), 'twoway');
fieldname = fieldnames(files);
for f=1:numel(fieldname)
    files.(fieldname{f})(~set_filter) = [];
end
files.LUT = LUT(set_filter,:);

% ------------ Prepare files for learning curve analysis ------------

if isfield(cfg, 'learning_curve') && cfg.learning_curve %% && cfg.decoding.scheme ~= 3
    chunk = files.chunk(files.cvclass); % only on the same class as CV works on!!!
    chunk_numbers = unique(chunk);
    
    % Cooncatenate files of each learning curve step
    fieldname(ismember(fieldname, 'set')) = []; % keep files.set separately
    for i=2:numel(chunk_numbers)
        for f=1:numel(fieldname)
            files.(fieldname{f}) = [files.(fieldname{f}); ...
                files.(fieldname{f})(ismember(chunk, chunk_numbers(1:i)))]; % the order of concat is essential!!
        end
        files.set = [files.set; repmat(i, sum(ismember(chunk, chunk_numbers(1:i))), 1)];
    end
%     unique(files.set)
    
    % Delete original set
    fieldname = [fieldname; {'set'}];
    for f=1:numel(fieldname)
        files.(fieldname{f})(files.set==1) = []; % delete original set
    end
    files.set = files.set - 1;
end

% Read ROI files if needed
if strcmp(cfg.analysis, 'roi')
    files.mask =  cellfun(@(x,y) fullfile(get_path('project'), ...
        get_folder(cfg.subject, 'r'), 'fMRI', 'ROI', x, 'coregistered', ...
        [y cfg.roi.extension]), cfg.roi.mask.atlas, cfg.roi.mask.fname, 'un', 0);
end