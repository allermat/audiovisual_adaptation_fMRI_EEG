function [predata, postdata] = load_test_data(subject, datatype)

% Preallocate data and define deafult data storage
[predata, postdata] = deal(cell(2, 2));
adapt = {'radapt' 'ladapt'; 'radapt' 'ladapt'};

% Load pre-test sessions
datafolder = fullfile(get_path('project'), get_folder(subject, 'r', datatype, 'behavioural data', 'pretest'));
if ~iscell(datafolder)
    datafolder = {datafolder};
end
nsessions = numel(datafolder);
if nsessions > 4
    error('maximum number of sessions exceeded')
end
for ses=1:nsessions
    s = sum(ismember(subject.(datatype).pretest.match(1:ses), subject.(datatype).pretest.match{ses})); % session within adaptation 
    a = find(cellfun(@(x) x(1), upper(adapt(1,:))) == subject.(datatype).pretest.match{ses}); % adaptation side
    predata{s,a} = loadfile(datafolder(ses), {strcat(subject.id, '_run*.mat')}, subject.(datatype).pretest.runid(ses), {'data'});
end

% Load post-test sessions
for a=1:2 % adaptation side
    if ~isempty(subject.(datatype).posttest.(adapt{1,a}).session)
        datafolder = fullfile(get_path('project'), get_folder(subject, 'r', datatype, 'behavioural data', 'posttest', adapt{1,a}));
        if ~iscell(datafolder)
            datafolder = {datafolder};
        end
        nsessions = numel(datafolder);
        if nsessions > 2
            error('maximum number of sessions exceeded')
        end
        for ses=1:nsessions
            postdata{ses,a} = loadfile(datafolder(ses), {strcat(subject.id, '_run*.mat')}, subject.(datatype).posttest.(adapt{ses,a}).runid(ses), {'data'});
        end
    end
end