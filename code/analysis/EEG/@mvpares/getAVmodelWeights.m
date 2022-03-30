function weights = getAVmodelWeights(obj,varargin)
% Method for accessing AV model weights
%
% USAGE:
%   weigths = getAVmodelWeights(obj)
%   weigths = getAVmodelWeights(obj,'Name',Value)
% INPUT:
%   Required:
%       obj (object): mvpares object
%   'Name'-Value arguments:
%       genTime (string): indicating the generalization time. Possible 
%           values: 'tr','tr_x_tr' for trainin time and training time by 
%           training time respectively. Default: 'tr'. 
%       smooth (logical): whether to smooth the time course of AV estimates
%       unit (string): the unit of the output, 'deg' - degreees, 
%           'rad' - radians, default: 'rad'
% OUTPUT:
%   weights (struct): weights as fields of the array. Each field
%       contains a time series of weights for a certain condition, or a
%       time x time matrix of weights depending on input settings. 

% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

% Parsing input
p = inputParser;
validGenTimes = {'tr','tr_x_tr'};
validUnits = {'deg','rad'};
addRequired(p,'obj');
addParameter(p,'genTime','tr',@(x)any(validatestring(x,validGenTimes)));
addParameter(p,'smooth',false,@(x)validateattributes(x,{'logical'},...
    {'scalar'}));
addParameter(p,'unit','rad',@(x)any(validatestring(x,validUnits)));
parse(p,obj,varargin{:});
obj = p.Results.obj;
genTime = p.Results.genTime;
smooth = p.Results.smooth;
unit = p.Results.unit;

if strcmp(obj.level,'group')
    % Group level AV weigths are stored in the dataset
    
    % Checking if AV model weigths are part of the dataset
    if any(~ismember({'gen_AVmodelWeights','gen_AVmodelWeights_s'},obj.who))
        weights = [];
        warning('mvpares:getAVmodelWeights:reqestedDataNotPresent',...
            'The dataset does not contain AV model weights. ');
        return;
    end
    
    if strcmp(genTime,'tr_x_tr')
        % Check if genTime is valid given the result dataset
        if isvector(obj.getGenTimePoints)
            weights = [];
            warning('mvpares:getAVmodelWeights:reqestedDataNotPresent',...
                'The dataset is not generalized time x time! ');
            return;
        end
        if smooth
            weights = [];
            warning('mvpares:getAVmodelWeights:invalidSettings',...
                ['Can''t smooth if time x time generalized estimates ',...
                'are required! ']);
            return;
        end
        weights = obj.data.gen_AVmodelWeights;
    elseif strcmp(genTime,'tr')
        if smooth
            weights = obj.data.gen_AVmodelWeights_s;
        else
            weights = obj.data.gen_AVmodelWeights;
            % Extracting the diagonal if just the training timepoints are 
            % needed and the data is time x time generalized
            if ~isvector(obj.getGenTimePoints)
                weights = structfun(@diag,weights,'UniformOutput',false);
            end
        end
    end
        
else
    
    % Subject level AV weigths are computed on-line from the AV estimates
    estimates = obj.getAVmodelEstimates('genTime',genTime,'smooth',smooth);
    
    if isempty(estimates)
        warning('mvpares:getAVmodelWeights:reqestedDataNotPresent',...
            'Couldn''t find AV model estimates in the dataset, returning.');
        weights = [];
        return;
    end
    
    % Separating A and V betas, matched by alphabetical order
    estimateNames = fieldnames(estimates);
    estimatesCell = struct2cell(estimates);
    
    aBetas = estimatesCell(~cellfun(@isempty,regexp(estimateNames,'A_.*_beta')));
    vBetas = estimatesCell(~cellfun(@isempty,regexp(estimateNames,'V_.*_beta')));
    aCI = estimatesCell(~cellfun(@isempty,regexp(estimateNames,'A_.*_ci')));
    vCI = estimatesCell(~cellfun(@isempty,regexp(estimateNames,'V_.*_ci')));
    
    aBetaNames = estimateNames(~cellfun(@isempty,regexp(estimateNames,'A_.*_beta')));
    condNames = strrep(strrep(aBetaNames,'A_',''),'_beta','');
    
    % Computing Wavs
    funb2w = @(v,a) atan2(v,a);
    wavs = cellfun(funb2w,vBetas,aBetas,'UniformOutput',false);
    aCIhigh = cellfun(@plus,aBetas,aCI,'UniformOutput',false);
    vCIhigh = cellfun(@plus,vBetas,vCI,'UniformOutput',false);
    wavsCIhigh = cellfun(funb2w,vCIhigh,aCIhigh,'UniformOutput',false);
    wavsCI = cellfun(@minus,wavsCIhigh,wavs,'UniformOutput',false);
    
    wavNames = strcat(condNames,repmat({'_wav'},size(condNames)));
    ciNames = strcat(condNames,repmat({'_ci'},size(condNames)));
    
    weights = cell2struct(cat(1,wavs,wavsCI),cat(2,wavNames,ciNames));
    
end

% Converting to degrees if necessary
if strcmp(unit,'deg')
    weights = structfun(@degrees,weights,'UniformOutput',false);
end

end

