function mvparesObj = addStats(mvparesObj,var,varargin)
% Method for computing various statistics
%
% USAGE:
%   mvparesObj = addStats(mvparesObj,var)
%   mvparesObj = addStats(mvparesObj,var,'Name',Value)
% INPUT:
%   Required:
%       mvparesObj (object): mvpares object
%       var (string): the variable for which the statistics are computed.
%           Possible values: 'avWeights'
%   'Name'-Value arguments:
%       pathIndivFiles (cell array): full path of individual files for
%           statistics on group level data (must be specified for group 
%           level data)
% OUTPUT:
%   mvparesObj (object): mvpares object with the required statistics added. 
 
% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

% Parsing input
p = inputParser;
validVars = {'genPerf','recalIndex'};
addRequired(p,'mvparesObj',@(x) isa(x,'mvpares') && x.isvalid);
addRequired(p,'var',@(x) any(validatestring(x,validVars)));
addParameter(p,'pathIndivFiles',{},@(x) iscellstr(x));
addParameter(p,'timeWin',{},@(x) iscell(x));
parse(p,mvparesObj,var,varargin{:});
mvparesObj = p.Results.mvparesObj;
var = p.Results.var;
pathIndivFiles = p.Results.pathIndivFiles;
timeWin = p.Results.timeWin;

% Assign empty array and return if the result dataset is not
% generalized.
if strcmp(mvparesObj.level,'subj')
    warning('mvpa:addStats:datasetLevelMismatch',...
        ['The dataset''s level is ''subj''. Subject level statistics are ',...
        'not yet implemented. Returning.']);
    return;
else
    if isempty(pathIndivFiles)
        warning('mvpa:addStats:missingInput',...
            'pathIndivFiles must be specified for group level statistics. Returning.');
        return;
    else
        if any(~cellfun(@exist,pathIndivFiles))
            warning('mvpa:addStats:invalidInput',...
                'Each path in pathIndivFiles must point to an existing file. Returning.');
            return;
        end
    end
end

info = mvparesObj.getInfo;
% This makes sure that both Windows and Unix paths are recognized,
% regardless of the operation system the code executes on. 
temp = regexp(info.sourceFiles,'.*[\\|/](.*)\.mat','tokens','once');
sourceFileNames = cat(1,temp{:});
temp = regexp(pathIndivFiles,'.*[\\|/](.*)\.mat','tokens','once');
indivFileNames = cat(1,temp{:});

% Checking if all required datasets are specified
if any(~ismember(sourceFileNames,indivFileNames))
    warning('mvpa:addStats:missingDataset',...
        'At least one of the required individual datasets is not found. Returning. ');
    return;
end

% Loading individual objects
objList = cellfun(@mvpares,pathIndivFiles,'UniformOutput',false);

switch var
    case 'genPerf'
        fieldName = 'stat_perfEstimates';
        
        indivData = cellfun(@getGenPerfEstimates,objList,...
                            repmat({'genTime'},size(objList)),...
                            repmat({'tr'},size(objList)));
        stats = mvpa.compGenPerfStats(indivData,mvparesObj.getTrTimePoints);
        
    case 'recalIndex'
        fieldName = 'stat_recalIndex';
        
        indivData = cellfun(@getRecalIndex,objList,...
                            repmat({'genTime'},size(objList)),...
                            repmat({'tr'},size(objList)));
        if isempty(timeWin)
            stats = mvpa.compRecalIndexStats(indivData, ...
                                             mvparesObj.getTrTimePoints);
        else
            stats = mvpa.compRecalIndexStats(indivData, ...
                                             mvparesObj.getTrTimePoints,timeWin);
        end
end

% Setting the dataset object to writable if it is not
if ~mvparesObj.writable, mvparesObj.setWritable(true); end
if ismember(fieldName,fieldnames(mvparesObj.data))
    warning('mvpa:addStats:overwriteField',...
            ['The field %s already exists in the mvpa result dataset, '...
             'it will be overwritten.'],fieldName);
end
mvparesObj.data.(fieldName) = orderfields(stats);
mvparesObj.setWritable(false);

end

 