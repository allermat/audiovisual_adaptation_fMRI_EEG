function LUT = mvpa_SPM_read(cfg)

% Spm folder
spmdir = fullfile(get_path('project'), get_folder(cfg.subject, 'r'), 'fMRI', ...
    '1st level', cfg.spmsubdir);

% Extract regressors from SPM
fprintf('extract regressors from SPM...');
load(fullfile(spmdir, 'SPM.mat'));

% Filter to hrf convolved aloc regressors
id = find(~cellfun(@isempty, regexp(SPM.xX.name, '.*aloc.*bf\(1\)$')));
betafname = arrayfun(@(x) sprintf('beta_%04d.nii', x), id, 'un', 0); % nii
regnames = SPM.xX.name(id);
LUT = table();
for n=1:numel(regnames)
    tmp = regexpi(regnames{n}, '^Sn\(([0-9]+)\)', 'tokens'); % get SPM session
    spm_session = cellfun(@str2double, tmp{:});
    tmp = regexpi(regnames{n}, '([a-z]+=-?\d+)', 'match'); % get factors and levels in <factor>=<level> format
    Tn = cell2table([spm_session cellfun(@(x) str2double(regexpi(x, '(-?\d+)', 'match')), tmp, 'un', 0)], ...
        'VariableNames', [{'spm_session'} cellfun(@(x) regexpi(x, '([a-z]+)', 'match'), tmp)]);
    LUT = [LUT; Tn];
end
% LUT.beta = fullfile(spmdir, betafname)';
LUT.betaid = id';
fprintf('done\n');
norigvars = size(LUT, 2);

% Additional operations on sessions, runs and blocks/periods
LUT(~ismember(LUT.session, cfg.subject.fMRI.session),:) = []; % remove sessions we do not want
if strfind(cfg.results.dir, 'multivar')
    isrowid = 0;
else
    isrowid = 1;
end
if strncmp(cfg.spmsubdir, 'run_merged', 10)
    tmp = table2array(unique(LUT(:,{'session' 'run'})));
    LUT(:,'overall_run') = rowfun(@(x,y) find(ismember(tmp, [x y], 'rows')), LUT, 'InputVariables', {'session' 'run'});
elseif strcmp(cfg.spmsubdir, 'merged')
    LUT.overall_run = LUT.spm_session;
    LUT = merge_session(LUT, isrowid);
else
    tmp = session_tuple(cfg.subject);
    LUT(:,'overall_run') = rowfun(@(x,y) find(cellfun(@(z) ismember([x y], z, 'rows'), tmp)), LUT, 'InputVariables', {'session' 'run'}); % define unique runs (taking into account merging!!!)
    LUT = merge_session(LUT, isrowid); % prepare data merging belonging to same runs
end
LUT(:,'blocktype') = rowfun(@(x,y) get_blocktype(cfg.subject, x, y), LUT, 'InputVariables', {'session' 'run'});

% Chunks used for cross-validation
nfolds = 4;
id = {ismember(LUT.blocktype, 'pretest'); ~ismember(LUT.blocktype, 'pretest')};
for i=1:size(id, 1)
    chunk = unique(LUT.overall_run(id{i,1}));
    chunk_cv = mod(1:numel(chunk), nfolds) + 1; % no randomization with enabling sessions most evenly distributed across folds
    idrow = arrayfun(@(x) chunk_cv(ismember(chunk, x)), LUT.overall_run(id{i,1}));
    LUT.chunk(id{i,1},1) = idrow;
end

% ----------------- Auxiliary functions ------------------

function blocktype = get_blocktype(subj, session, run)
%   Gets blocktype (pretest, posttest, l-adaptation, r-adaptation) for a 
% subject's given session and run

blocktype = '';
field = {{'pretest'} {'posttest' 'radapt'} {'posttest' 'ladapt'}};
for i=1:numel(field)
    S = getfield(subj.fMRI, {1}, field{i}{:});
    runid = S.runid(S.session == session);
    if ~isempty(runid) && ismember(run, runid{:})
        blocktype = {strjoin(field{i}, '-')};
        return
    end
end
if isempty(blocktype)
   error('session %d and run %d not found in subject data', session, run); 
end


function tuple = session_tuple(subj)

% Create session-run tuples
tuple = [];
for s=1:numel(subj.fMRI.pretest.session)
    for r=1:numel(subj.fMRI.pretest.runid{s})
        tuple{end+1} = [subj.fMRI.pretest.session(s) subj.fMRI.pretest.runid{s}(r)];
    end
end
adapt = {'radapt' 'ladapt'};
for a=1:numel(adapt)
    for s=1:numel(subj.fMRI.posttest.(adapt{a}).session)
        for r=1:numel(subj.fMRI.posttest.(adapt{a}).runid{s})
            if ismember(r, subj.fMRI.posttest.(adapt{a}).mergerun{s})
                tuple{end} = [tuple{end}; [subj.fMRI.posttest.(adapt{a}).session(s) subj.fMRI.posttest.(adapt{a}).runid{s}(r)]];
            else
                tuple{end+1} = [subj.fMRI.posttest.(adapt{a}).session(s) subj.fMRI.posttest.(adapt{a}).runid{s}(r)];
            end
        end
    end
end

% Sort tuple
tmp = cellfun(@(x) x(1,:), tuple, 'un', 0);
[B, I] = sortrows(cat(1, tmp{:}), 1);
tuple = tuple(I);


function LUT_nodup = merge_session(LUT, isrowid)

tuple = table2array(unique(LUT(:,{'overall_run' 'resp' 'aloc'})));
LUT_nodup = table();
betaid = [];
if isrowid
    rowid = [];
end
for i=1:size(tuple, 1)
   id = find((LUT.overall_run==tuple(i,1) & LUT.resp==tuple(i,2) & LUT.aloc==tuple(i,3)));
   tmp = LUT(id,:);
   LUT_nodup = [LUT_nodup; tmp(1,~ismember(tmp.Properties.VariableNames, 'betaid'))]; 
   betaid(end+1,1:2) = tmp.betaid;
   if isrowid
       rowid(end+1,1:2) = id;
   end
end
LUT_nodup.betaid = betaid;
if isrowid
    LUT_nodup.rowid = rowid;
end