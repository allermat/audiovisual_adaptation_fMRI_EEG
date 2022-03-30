function out = run_behav_all(modality)

% Check modality
if ~ismember(modality, {'psychophysics' 'fMRI' 'EEG'})
    error('Input argument should be one of the following modalities: ''psychophysics'' ''fMRI'' ''EEG''');
end

% Load subject info
subjects = subject_info;
if strcmp(modality, 'psychophysics')
    nsubjects = 15;
else
    nsubjects = 5;
end

% AV-adaptation (radapt) and VA-adaptation (ladapt)
adapt = {'radapt' 'ladapt'};
    
for ss=1:nsubjects
    % Load pre- and post-test data
    [out(ss).pre.data, out(ss).post.data] = load_test_data(subjects(ss), modality);
    
    % Calculate number of positive and out of number responses
    [out(ss).pre.NumPos, out(ss).pre.OutOfNum] = get_numpos_outofnum(out(ss).pre.data);
    [out(ss).post.NumPos, out(ss).post.OutOfNum] = get_numpos_outofnum(out(ss).post.data);
    
    % Calculate adaptation phase hit rate, false alarm rate and d-prime
    for a=1:2
        funcfolder = fullfile(get_path('project'), get_folder(subjects(ss), ...
            'r', modality, 'behavioural data', 'posttest', adapt{a}));
        if ~iscell(funcfolder)
            funcfolder = {funcfolder};
        end
        for s=1:numel(funcfolder)
            varargin = {'exp', modality, 'subjid', subjects(ss).id, ...
                'session', subjects(ss).(modality).posttest.(adapt{a}).session(s), ...
                'runid', subjects(ss).(modality).posttest.(adapt{a}).runid{s}, ...
                'plot', 0, 'print', 0, 'blockid', repmat({[1:2:10]}, 1, ...
                (numel(subjects(ss).(modality).posttest.(adapt{a}).runid{s})))};
            if ss < 6
                varargin(end+1:end+2) = {'session_folder_name', subjects(ss).session_folder_name};
            end
            [hit{a,s}, false_alarm{a,s}] = subject_analysis_catch_trial('offline', varargin{:});
        end
        out(ss).adapt.dprime.(adapt{a}) = calc_dprime(mean(cat(2, hit{a,:})), ...
            mean(cat(2, false_alarm{a,:})));
        out(ss).adapt.hit.(adapt{a}) = mean(cat(2, hit{a,:})) * 100;
        out(ss).adapt.false_alarm.(adapt{a}) = mean(cat(2, false_alarm{a,:})) * 100;
    end
    clear hit false_alarm
    
    % Calculate pre-adaptation phase hit rate, false alarm rate and d-prime
    funcfolder = fullfile(get_path('project'), get_folder(subjects(ss), ...
        'r', modality, 'behavioural data', 'pretest'));
    if ~iscell(funcfolder)
        funcfolder = {funcfolder};
    end
    for s=1:numel(funcfolder)
        varargin = {'exp', modality, 'subjid', subjects(ss).id, ...
            'session', subjects(ss).(modality).pretest.session(s), ...
            'runid', subjects(ss).(modality).pretest.runid{s}, 'plot', 0, 'print', 0};
        if ss < 6
            varargin(end+1:end+2) = {'session_folder_name', subjects(ss).session_folder_name};
        end
        [hit{s}, false_alarm{s}] = subject_analysis_catch_trial('offline', varargin{:});
    end
    out(ss).pre.dprime = calc_dprime(mean(cat(2, hit{:})), mean(cat(2, false_alarm{:})));
    out(ss).pre.hit = mean(cat(2, hit{:})) * 100;
    out(ss).pre.false_alarm = mean(cat(2, false_alarm{:})) * 100;
    clear hit false_alarm
    
    % Calculate post-adaptation phase hit rate, false alarm rate and d-prime
    for a=1:2
        funcfolder = fullfile(get_path('project'), get_folder(subjects(ss), ...
            'r', modality, 'behavioural data', 'posttest', adapt{a}));
        if ~iscell(funcfolder)
            funcfolder = {funcfolder};
        end
        for s=1:numel(funcfolder)
            varargin = {'exp', modality, 'subjid', subjects(ss).id, ...
                'session', subjects(ss).(modality).posttest.(adapt{a}).session(s), ...
                'runid', subjects(ss).(modality).posttest.(adapt{a}).runid{s}, ...
                'plot', 0, 'print', 0, 'blockid', repmat({[2:2:10]}, 1, ...
                    (numel(subjects(ss).(modality).posttest.(adapt{a}).runid{s})))};
            if ss < 6
                varargin(end+1:end+2) = {'session_folder_name', subjects(ss).session_folder_name};
            end
            [hit{a,s}, false_alarm{a,s}] = subject_analysis_catch_trial('offline', varargin{:});
        end
        out(ss).post.dprime.(adapt{a}) = calc_dprime(mean(cat(2, hit{a,:})), ...
            mean(cat(2, false_alarm{a,:})));
        out(ss).post.hit.(adapt{a}) = mean(cat(2, hit{a,:})) * 100;
        out(ss).post.false_alarm.(adapt{a}) = mean(cat(2, false_alarm{a,:})) * 100;
    end
    clear hit false_alarm
end