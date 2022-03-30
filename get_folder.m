function datafolder = get_folder(subject, permission, varargin)
%   get_folder builds up the relative path to any data folder
%
% NOTES:
%   1. folder is automatically created with writing permission (w), input is
%   required with reading permission (r) if the folder does not exist
%   2. only existing session folder can be requested with read permission
%   3. existing session folder can be requested with write permission only
%       on the day of creation

% Main data folder
datafolder = 'data';

% Group or subject data folder
if isstruct(subject) && isfield(subject, 'id') % subject
    if isfield(subject, 'pilot') && ~isempty(subject.pilot)
        datafolder = fullfile(datafolder, sprintf('pilot_%d', subject.pilot), subject.id);
    else
        datafolder = fullfile(datafolder, subject.id);
    end
    folders = getfname(fullfile(get_path('project'), datafolder), '*');
    
    imaging = {'fMRI' 'EEG'};
    if any(ismember(imaging, folders))
        imaging = imaging(ismember(imaging, folders));
        for i=1:numel(imaging)
            folders{ismember(folders, imaging(i))} = fullfile(imaging{i}, getfname(fullfile(get_path('project'), datafolder, imaging{i}), '*'));
            folders = cat(1, folders{:});
        end
    end
    % List any files and folders in data folders (of any task)
    if ~isempty(folders)
        if isfield(subject, 'session_folder_name') && ...
                ~isempty(subject.session_folder_name)
            session = subject.session_folder_name;
        else
            for f=1:numel(folders)
                fname{f} = getfname(fullfile(get_path('project'), datafolder, folders{f}), '*');
            end
            fname = unique(cat(1, fname{:}));
            % Find sessions based on format
            expression = {'(\d{4})_(\d{2})_(\d{2})'}; % session format e.g. 2016_11_12
            session = regexp(fname, expression);
            session = fname(~cellfun(@isempty, session));
        end
    else
        session = [];
    end
    datafolder = fullfile(datafolder, varargin{:});
    [strpath, name, ext] = fileparts(datafolder);
    while ~isdir(fullfile(get_path('project'), datafolder)) && strcmp(permission, 'r') && isempty(str2num(name)) % find data folder recursively
        [datafolder, name, ext] = fileparts(datafolder);
    end
    datafolder = {datafolder}; % change to cell mode...:)
    fieldname = varargin;
    i = 1;
    S = subject;
    if isfield(S, 'session')
        ses = S.session;
    else
        while i <= numel(fieldname) % find session info recursively
            if isfield(S, fieldname{i})
                S = getfield(subject, {1}, fieldname{1:i});
                i = i + 1;
            else
                fieldname(i) = [];
            end
        end
        if isfield(S, 'session')
            ses = S.session;
        end
    end
    if isempty(session) && strcmp(permission, 'w')
        datafolder = fullfile(datafolder, datestr(now, 'yyyy_mm_dd'));
    elseif exist('ses', 'var')
        if any(ses > length(session)) && strcmp(permission, 'r')
            error('session folder does not exist');
        elseif any(ses <= length(session)) % session already exists
            if strcmp(permission, 'w') && ~strcmp(session{ses}, datestr(now, 'yyyy_mm_dd')) % session date is different!!
                warning('given session folder already exists with another date');
            end
            datafolder = fullfile(datafolder, session(ses));
        else
            datafolder = fullfile(datafolder, datestr(now, 'yyyy_mm_dd'));
        end
        if isempty(datafolder)
           error('session folder not found'); 
        end
        if isfield(S, 'runid') && ~ismember('psychophysics', varargin)
            for s=1:numel(ses)
                datafolder{s} = fullfile(datafolder{s}, arrayfun(@(x) getfname(fullfile(get_path('project'), datafolder{s}), sprintf('*run%02d*', x)), S.runid{s}')); % ses(s)
            end
            datafolder = cat(1, datafolder{:});
            id = cellfun(@(x) isdir(fullfile(get_path('project'), x)), datafolder); % true folders
            datafolder = [datafolder(id); unique(cellfun(@fileparts, datafolder(~id), 'un', 0))]; % % filter for true folders
        end
    end
elseif strcmp(subject, 'group') % group
    datafolder = fullfile(datafolder, {'group'}, varargin{:});
end

% Create folder if needed
for d=1:numel(datafolder)
    if ~isdir(fullfile(get_path('project'), datafolder{d}))
        switch permission
            case 'r'
                disp(fullfile(get_path('project'), datafolder{d}));
                response = input('data folder does not exist! do you want to create it? (y/n) ', 's');
                switch response
                    case 'y'; mkdir(fullfile(get_path('project'), datafolder{d}));
%                     case 'n'; error('program terminated');
                end
            case 'w'
                mkdir(fullfile(get_path('project'), datafolder{d}));
        end
    end
end
if numel(datafolder) == 1
    datafolder = datafolder{:}; % convert to char array
end
