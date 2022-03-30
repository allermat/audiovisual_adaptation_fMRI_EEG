function mvparesObj = addNeurometricFunctions(mvparesObj)
% Method for computing neurometric functions
% 
% USAGE:
%   mvparesObj = addGenPerfEstimates(mvparesObj)
% INPUT:
%   mvparesObj (object): mvpares object
% OUTPUT:
%   mvparesObj (object): mvpares object with genaralization performance 
%       estimates

% Copyright(C) 2017, Mate Aller
% allermat@gmail.com

% Parsing input
p = inputParser;
addRequired(p,'mvparesObj',@(x) isa(x,'mvpares') && x.isvalid);
parse(p,mvparesObj);
mvparesObj = p.Results.mvparesObj;

sGenTime = mvparesObj.getSizeGenTime;
if sGenTime(2) > 1
    genTimeStr = 'tr_x_tr';
else
    genTimeStr = 'tr';
end

% Extract necessary data
predLabels = mvparesObj.getPredLabels('genTime',genTimeStr);
infoGenExamples = mvparesObj.getInfoGenExamples;
if isempty(predLabels) || isempty(infoGenExamples)
    warning('mvpa:addNeurometricFunctions:requestedDataNotPresent',...
        ['Couldn''t extract necessary data, neurometric functions ',...
        'can''t be estimated, returning.']);
    return;
end

NMF = mvpa.fitNeurometricFunctions(predLabels,infoGenExamples,genTimeStr);

% Setting the dataset object to writable if it is not
if ~mvparesObj.writable, mvparesObj.setWritable(true); end

fieldName = 'gen_neuroMetrFun';
if ismember(fieldName,fieldnames(mvparesObj.data))
    warning('mvpa:addNeurometricFunctions:overwriteField',...
        ['The field ''%s'' already exists in the ',...
        'mvpa result dataset, it will be overwritten. '],fieldName);
end
mvparesObj.data.(fieldName) = NMF;
mvparesObj.setWritable(false);

end