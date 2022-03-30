function mvparesObj = addRecalIndex(mvparesObj)
% Method for computing the recalibration index
% 
% DETAILS: 
%   RI = mean(A_pred_R) - mean(A_pred_L)
% 
% USAGE:
%   mvparesObj = addRecalIndex(mvparesObj)
% INPUT:
%   mvparesObj (object): mvpares object
% OUTPUT:
%   mvparesObj (object): mvpares object with recalibration index

% Copyright(C) 2016, Mate Aller
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

predLabels = mvparesObj.getPredLabels('genTime',genTimeStr);
predLabelInfo = mvparesObj.getInfoGenExamples;
if isempty(predLabels) || isempty(predLabelInfo)
    warning('mvpa:addRecalIndex:requestedDataNotPresent',...
        ['Couldn''t extract predicted labels, performance estimates ',...
        'can''t be estimated, returning.']);
    return;
end

[temp.RI,temp.RIbin] = mvpa.compRecalIndex(predLabelInfo,predLabels);

% Setting the dataset object to writable if it is not
if ~mvparesObj.writable, mvparesObj.setWritable(true); end

fieldName = 'gen_recalIndex';
if ismember(fieldName,fieldnames(mvparesObj.data))
    warning('mvpa:addRecalIndex:overwriteField',...
        ['The field ''%s'' already exists in the ',...
        'mvpa result dataset, it will be overwritten. '],fieldName);
end
mvparesObj.data.(fieldName) = temp;
mvparesObj.setWritable(false);

end