function [hit_rate, false_alarm_rate] = subject_analysis_catch_trial(analmode, varargin)
% subject_analysis_catch_trial(analmode, varargin)
%   subject level analysis either for analysis on the fly (see online mode,
%   this way the function is called by present_stimuli.m after a completed
%   block) or offline data check
%
% Input:
%   analmode = online or offline
%
% Examples:
% 1. online
%   subject_analysis_catch_trial('online', data)
%
% 2. offline
%   subject_analysis_catch_trial('offline', 'subjid', '334', 'session', 2, 'runid', 1, 'blockid', {1:2})
%
% Notes:
% 1. In online mode, data should be given in input without any other input
%   arguments. Blockid is not relevant, because we always want to analyse
%   all the data that is given.
%
% 2. In offline mode, data should be loaded, therefore the many input
%   arguments. Some default values are also defined.
% 
% 3. At the moment, plotting a figure is not supported yet.
%
% 16/11/2016 AM

if nargin == 0
   error('input arguments missing, see help') 
end

% Load data and handle inputs based on analysis mode
switch analmode
    case 'online'        
        % Assign data
        if nargin ~= 2
            error('input variables should be given as in the example of the help')
        end
        if strcmp(inputname(2), 'data')
            data = varargin{1};
        else
            error('data should be given in online mode, see help')
        end
        
        % Instead of input args we define values
        S = struct('runid', unique(data.run), 'figure', 0, 'print', 0);
        
    case 'offline'
        % Assign possible input arguments (overwriting default values)
        if ~mod(nargin, 2)
            error('variable input arguments should come in pairs')
        end
        S = struct('exp', 'psychophysics', 'figure', 0, 'print', 1);
        for i=1:2:numel(varargin)
            S.(varargin{i}) = varargin{i+1};
        end
        
        % Error handling
        if any(~isfield(S, {'subjid' 'session' 'runid'}))
            error('subjid, session and runid should be given, see help');
        end
        
        if isfield(S, 'blockid')
            flag = ['b' sprintf('_%d', S.blockid{:})];
        else
            flag = 'b_All';
        end
      
        % Load data
        for d=1:numel(S.runid)
            if isfield(S, 'pilot')
                subject = struct('id', S.subjid, 'session', S.session, 'pilot', S.pilot);
            else
                subject = struct('id', S.subjid, 'session', S.session);
            end
            if isfield(S, 'session_folder_name')
                subject.session_folder_name = S.session_folder_name;
            end
            datafolder = fullfile(get_path('project'), get_folder(subject, 'r', S.exp, 'behavioural data'));
            try fname = arrayfun(@(x) getfname(datafolder, sprintf('*run%02d*', x)), S.runid(d)); catch; error('file not found'); end
            data{d} = load_truncated_data(fullfile(datafolder, fname{1}), 20, 'data');
%             data{d} = loadfile({datafolder}, {strcat(subject.id, '_run*.mat')}, {S.runid(d)}, {'data'});
            if isfield(S, 'blockid') % only blocks of data we are interested in
                data{d} = data{d}(ismember(data{d}.block, S.blockid{d}),:);
            end
        end
        data = cat(1, data{:});
  
end

% -------------- Calculate hit and false alarm rate ----------------

if sum(data.catch_trial) == 0
   fprintf('no catch trials in the given data\n');
   [hit_rate, false_alarm_rate] = deal(NaN);
   return;
end

% Define overall block id based on runs and blocks
try
    id = unique(double(data(:,{'run' 'block'})), 'rows');
catch
    id = [1 1];
    data.run = ones(size(data, 1), 1);
    data.block = ones(size(data, 1), 1);
end
ntrials = size(data, 1);
data.id = zeros(ntrials, 1);
for i=1:ntrials
    data.id(i) = find(data.run(i) == id(:,1) & data.block(i) == id(:,2));
end

block = unique(data.id);
nblocks = numel(block);

for b=1:nblocks
    blockdata = data(data.id==block(b),:);
    catch_trial = blockdata.catch_trial;    
    hit_rate(b) = sum(catch_trial & ~isnan(blockdata.resp)) / sum(catch_trial);
    false_alarm_rate(b) = sum(~catch_trial & ~isnan(blockdata.resp)) / sum(catch_trial);
    dprime_all(b) = calc_dprime(hit_rate(b), false_alarm_rate(b));
    if ~isempty(cell2mat(strfind(unique(blockdata.blocktype), 'test')))
        tasktype{b} = 'alocalization';
    elseif ~isempty(cell2mat(strfind(unique(blockdata.blocktype), 'adapt')))
        tasktype{b} = 'vdetection';
    end
end