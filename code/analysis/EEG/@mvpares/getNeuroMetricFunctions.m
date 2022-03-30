function nmf = getNeuroMetricFunctions(obj,varargin)
% Method for accessing neurometric functions
%
% USAGE:
%   nmf = getNeurometricFunctions(obj)
%   nmf = getNeurometricFunctions(obj,'Name',Value)
% INPUT:
%   Required:
%       obj (object): mvpares object
%   'Name'-Value arguments:
%       genTime (string): indicating the generalization time. Possible 
%           values: 'tr','tr_x_tr' for trainin time and training time by 
%           training time respectively. Default: 'tr'. 
%       inference (string): group level inference 'rfx' for random
%           effects analysis, 'ffx' for fixed effects analysis
% OUTPUT:
%   nmf (struct): neurometric function parameters as fields of the
%       structure. 
%       Fields: 
%           PF: function handle for the fitted function
%           stimLevels: stimulus location levels
%           *_pv: array of fitted parameter values for the
%               given condition '*' (threshold,slope,guess rate,lapse rate)
%           *_pctr: array of % decoded right for the given
%               condition '*'
%       

% Copyright(C) 2017, Mate Aller
% allermat@gmail.com

% Parsing input
p = inputParser;
validGenTimes = {'tr','tr_x_tr'};
validInference = {'ffx','rfx'};
addRequired(p,'obj');
addParameter(p,'genTime','tr',@(x)any(validatestring(x,validGenTimes)));
addParameter(p,'inference','tr',@(x)any(validatestring(x,validInference)));
parse(p,obj,varargin{:});
obj = p.Results.obj;
genTime = p.Results.genTime;
inference = p.Results.inference;

% Assign empty array and return if the result dataset is not
% generalized.
if strcmp(obj.state,'trained')
    nmf = [];
    warning('mvpares:getNeurometricFunctions:datasetStateMismatch',...
        ['The dataset''s state is ''traned'', so it does ',...
        'not contain generalization data']);
    return;
end

% Assign empty array and return if the result dataset does not contain
% neurometric functions
if ~ismember('gen_neuroMetrFun',obj.who)
    nmf = [];
    warning('mvpares:getNeurometricFunctions:reqestedDataNotPresent',...
        'The dataset does not contain AV model estimates. ');
    return;
end

if strcmp(inference,'ffx')
    nmf = obj.data.gen_neuroMetrFun_ffx;
else
    nmf = obj.data.gen_neuroMetrFun;
end

if strcmp(genTime,'tr_x_tr')
    % Check if genTime is valid given the result dataset
    if isvector(obj.getGenTimePoints)
        nmf = [];
        warning('mvpares:getNeurometricFunctions:reqestedDataNotPresent',...
                'The dataset is not generalized time x time! ');
        return;
    end
elseif strcmp(genTime,'tr')
    % Extracting the diagonal if just the training timepoints are needed
    % and the data is time x time generalized
    if ~isvector(obj.getGenTimePoints)
        temp = fieldnames(nmf);
        fieldsToRm = temp(cellfun(@isempty,regexp(temp,'(_pv|_pctr)$','once')));
        temp = rmfield(nmf,fieldsToRm);
        if strcmp(obj.level,'group') && strcmp(inference,'rfx')
            fun = @(x) reshape(x(:,repmat(logical(eye(size(x,2))),1,1,size(x,4))),...
                               size(x,1),size(x,2),size(x,4));
        else
            fun = @(x) x(:,logical(eye(size(x,2))));
        end
        temp = structfun(fun,temp,'UniformOutput',false);
        for i = 1:numel(fieldsToRm)
            temp.(fieldsToRm{i}) = nmf.(fieldsToRm{i});
        end
        nmf = temp;
    end
end


end