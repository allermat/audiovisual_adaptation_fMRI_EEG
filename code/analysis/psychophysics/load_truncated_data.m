function varargout = load_truncated_data(fname, ntrials, varargin)

% Load data
S = load(fname, varargin{:}); % load specific fields if provided

% Find block starts
if ismember('blockstart', S.data.Properties.VarNames)
    blockstart = [find(S.data.blockstart); size(S.data, 1)+1];
else
    warning('%s could not be truncated due to old data format', fname);
end

% Truncate data if needed
if diff(blockstart(end-1:end)) < ntrials
    id = 1:blockstart(end-1)-1;
    warning('%s data truncated with %d trials', fname, diff(blockstart(end-1:end)));
    if isempty(varargin)
        varargout{1} = S.data(id,:); % we assume that data should be truncated
    else
        for i=1:numel(varargin) % variables in varargin will be truncated
            switch varargin{i}
                case 'data'
                    varargout{i} = S.data(id,:);
                case 'stim'
                    S.stim.visual.loc = S.stim.visual.loc(:,:,id);
                    S.stim.audio.loc = S.stim.audio.loc(:,:,id);
                    varargout{i} = S.stim;
                case 'time'
                    field = {'predvonset' 'vonset' 'voffset' 'aonset' 'aoffset' 'predaonset'};
                    for f=1:numel(field)
                        S.time.(field{f}) = S.time.(field{f})(id);
                    end
                    varargout{i} = S.time;
            end
        end
    end
else
    fieldname = fieldnames(S);
    for i=1:numel(fieldname)
       varargout{i} = S.(fieldname{i});
    end
end