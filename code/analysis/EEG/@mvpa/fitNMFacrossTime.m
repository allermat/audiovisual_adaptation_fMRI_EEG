function varargout = fitNMFacrossTime(PF,stimLevels,numPos,outOfNum,varargin)
% Method to fit neurometric functions across training and
% generalization samples.

% Parsing input
p = inputParser;
addRequired(p,'PF',@(x) isa(x,'function_handle'));
addRequired(p,'stimLevels',@(x) validateattributes(x,{'numeric'}, ...
                                                  {'column'}));
addRequired(p,'numPos',@(x) validateattributes(x,{'numeric'}, ...
                                               {'nrows',size(stimLevels,1)}));
addRequired(p,'outOfNum',@(x) validateattributes(x,{'numeric'}, ...
                                               {'nrows',size(stimLevels,1)}));
addOptional(p,'opt',struct(),@(x) isstruct(x));
parse(p,PF,stimLevels,numPos,outOfNum,varargin{:});
PF = p.Results.PF;
stimLevels = p.Results.stimLevels;
numPos = p.Results.numPos;
outOfNum = p.Results.outOfNum;
opt = p.Results.opt;

[~,nTrTimePoints,nGenTimePointsPerTrTimePoint] = size(numPos);

pv = NaN(4,nTrTimePoints,nGenTimePointsPerTrTimePoint);
pctr = NaN(size(stimLevels,1),nTrTimePoints,nGenTimePointsPerTrTimePoint);

parfor i = 1:nTrTimePoints
    for j = 1:nGenTimePointsPerTrTimePoint
        actNumPos = numPos(:,i,j);
        paramsValues = fit_PF2(PF,stimLevels,actNumPos,outOfNum,opt);
        pv(:,i,j) = paramsValues;
        pctr(:,i,j) = actNumPos./outOfNum;
    end
end

varargout{1} = pv;
varargout{2} = pctr;

end