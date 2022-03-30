function corrCoeffs = getAVmodelCorrelations(obj,var,varargin)
% Method for accessing correlation coefficients from AV model correlations
%
% USAGE:
%   corrCoeffs = getAVmodelCorrelations(obj,var)
%   corrCoeffs = getAVmodelCorrelations(obj,var,'Name',Value)
% INPUT:
%   Required:
%       obj (object): mvpares object
%       var (string): the variable whos AV weights were correlated with 
%           the AV weights of the mvpares object. Possible values:
%           'acrossTime', 'behav', 'fmri'
%   'Name'-Value arguments:
%       genTime (string): indicating the generalization time. Possible 
%           values: 'tr','tr_x_tr' for trainin time and training time by 
%           training time respectively. Default: 'tr'. 
%       smooth (logical): whether to smooth the time courses
% OUTPUT:
%   corrCoeffs (struct): correlation coefficients as fields of the array.
%       Each field contains a time series of estimates of a certain 
%       condition, or a time x time matrix of estimates depending on input 
%       settings. 

% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

% Parsing input
p = inputParser;
validVars = {'all','acrossTime','behav','fmri'};
validGenTimes = {'tr','tr_x_tr'};
addRequired(p,'obj');
addRequired(p,'var',@(x) any(validatestring(x,validVars)));
addParameter(p,'genTime','tr',@(x)any(validatestring(x,validGenTimes)));
addParameter(p,'smooth',false,@(x)validateattributes(x,{'logical'},...
    {'scalar'}));
parse(p,obj,var,varargin{:});
obj = p.Results.obj;
genTime = p.Results.genTime;
var = p.Results.var;
smooth = p.Results.smooth;

% Assign empty array and return if the result dataset is not
% generalized.
if strcmp(obj.state,'trained')
    corrCoeffs = [];
    warning('mvpares:getAVmodelCorrelations:datasetStateMismatch',...
        ['The dataset''s state is ''traned'', so it does ',...
        'not contain generalization data']);
    return;
end

% Assign empty array and return if the result dataset does not contain
% AV model estimates
if any(~ismember({'gen_AVmodelCorrelations','gen_AVmodelCorrelations_s'},obj.who))
    corrCoeffs = [];
    warning('mvpares:getAVmodelCorrelations:reqestedDataNotPresent',...
        'The dataset does not contain AV model correlations. ');
    return;
end

% Loading estimates
if smooth
    corrCoeffs = obj.data.gen_AVmodelCorrelations_s;
else
    corrCoeffs = obj.data.gen_AVmodelCorrelations;
end
% Return all available data if required. Be aware that in this case the
% genTime setting might not be valid. Use with caution. 
if strcmp(var,'all'), return; end

% Getting rid of not required fields
allFields = fieldnames(corrCoeffs);
fieldsNotRequired = allFields(cellfun(@isempty,regexp(allFields,['.*_',var],'once')));
corrCoeffs = rmfield(corrCoeffs,fieldsNotRequired);

if strcmp(var,'acrossTime')
    if strcmp(genTime,'tr')
        % Extracting the diagonal if just the training timepoints are needed
        % and the data is time x time generalized
        corrCoeffs = structfun(@diag,corrCoeffs,'UniformOutput',false);
    end
else
    if strcmp(genTime,'tr_x_tr')
        % Check if genTime is valid given the result dataset
        if isvector(obj.getGenTimePoints)
            corrCoeffs = [];
            warning('mvpares:getAVmodelCorrelations:reqestedDataNotPresent',...
                'The dataset is not generalized time x time! ');
            return;
        elseif smooth
            corrCoeffs = [];
            warning('mvpares:getAVmodelCorrelations:reqestedDataNotPresent',...
                ['This correlation is not available smoothed in ',...
                'time x time generalized form. ']);
            return;
        end
    elseif strcmp(genTime,'tr')
        % Extracting the diagonal if just the training timepoints are needed
        % and the data is time x time generalized
        if ~isvector(obj.getGenTimePoints) && ~smooth
            corrCoeffs = structfun(@diag,corrCoeffs,'UniformOutput',false);
        end
    end
end

end