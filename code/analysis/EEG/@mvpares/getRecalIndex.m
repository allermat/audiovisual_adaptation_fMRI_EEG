function RI = getRecalIndex(obj,varargin)
% Method for accessing recalibration index
%
% USAGE:
%   RI = getRecalIndex(obj)
%   RI = getRecalIndex(obj,'Name',Value)
% INPUT:
%   Required:
%       obj (object): mvpares object
%   'Name'-Value arguments:
%       genTime (string): indicating the generalization time. Possible 
%           values: 'tr','tr_x_tr' for trainin time and training time by 
%           training time respectively. Default: 'tr'. 
%       smooth (logical): whether to smooth the time course of RI
% OUTPUT:
%   RI (struct): recalibration index structure. Fields: RI, RI_sem

% Copyright(C) 2017, Mate Aller
% allermat@gmail.com

% Parsing input
p = inputParser;
validGenTimes = {'tr','tr_x_tr'};
validInference = {'rfx','ffx'};
addRequired(p,'obj');
addParameter(p,'genTime','tr',@(x)any(validatestring(x,validGenTimes)));
addParameter(p,'inference','rfx',@(x) any(validatestring(x,validInference)));
addParameter(p,'smooth',false,@(x)validateattributes(x,{'logical'},...
    {'scalar'}));
parse(p,obj,varargin{:});
obj = p.Results.obj;
genTime = p.Results.genTime;
inference = p.Results.inference;
smooth = p.Results.smooth;

% Assign empty array and return if the result dataset is not
% generalized.
if strcmp(obj.state,'trained')
    RI = [];
    warning('mvpares:getRecalIndex:datasetStateMismatch',...
        ['The dataset''s state is ''traned'', so it does ',...
        'not contain generalization data']);
    return;
end

% Assign empty array and return if the result dataset does not contain
% recalibration index
if strcmp(inference,'rfx')
    if ~ismember('gen_recalIndex',obj.who)
        RI = [];
        warning('mvpares:getRecalIndex:reqestedDataNotPresent',...
                'The dataset does not contain the recalibration index. ');
        return;
    else
        % Loading estimates
        RI = obj.data.gen_recalIndex;
    end
else
    if ~ismember('gen_recalIndex_ffx',obj.who)
        RI = [];
        warning('mvpares:getRecalIndex:reqestedDataNotPresent',...
                'The dataset does not contain fixed effects recalibration index. ');
        return;
    else
        % Loading estimates
        RI = obj.data.gen_recalIndex_ffx;
    end 
end

if strcmp(genTime,'tr_x_tr')
    % Check if genTime is valid given the result dataset
    if isvector(obj.getGenTimePoints)
        RI = [];
        warning('mvpares:getRecalIndex:reqestedDataNotPresent',...
            'The dataset is not generalized time x time! ');
        return;
    elseif smooth
        RI = [];
        warning('mvpares:getRecalIndex:invalidSettings',...
            ['Can''t smooth if time x time generalized estimates ',...
            'are required! ']);
        return;
    end
elseif strcmp(genTime,'tr')
    % Extracting the diagonal if just the training timepoints are needed
    % and the data is time x time generalized
    if ~isvector(obj.getGenTimePoints)
        RI = structfun(@diag,RI,'UniformOutput',false);
    end
    
    if smooth

        % Smoothing the time courses moving average window of 20 ms
        fs = obj.getFsample;
        mavgWinSec = 0.02;
        mavgWinSample = round(mavgWinSec*fs);
        % Kernel for smoothing trials with the given time wintow
        kern = ones(mavgWinSample,1)./mavgWinSample;
        % Moving average function
        funMavg = @(x) conv(x,kern,'same');
        % Smoothing
        RI = structfun(funMavg,RI,'UniformOutput',false);

    end 
end

end