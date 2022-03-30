function trialInfoOut = buildtrialinfo(trialInfo,behavData,fileSpec)
% Builds trialInfo from the available behavioural information
% 
% DETAILS:
%   Converts the trialinfo field provided by FieldTrip to a
%   cellarray of structures, which makes the various pieces of
%   information more intuitive to access.
%       
% INPUT:
%   trialInfo (matrix): the trialinfo field of the FieltTrip data
%       structure
%   behavData: sturcture array of behavioural data (each array
%       element contains data from one run)
%   fileSpec: structure with information about the EEG file
% 
% OUTPUT: 
%   trialInfoOut (cell): cell array of structures as containing
%       trial informaiton

% Parsing input
p = inputParser;

% Input checking functions
checkBehavData = @(x) isstruct(x) && all(isfield(x,{'run', ...
                    'data'}));
checkFileSpec = @(x) isstruct(x) && all(isfield(x,{'name','session',...
                   'nRunsInFile','runid','hand','runsToExclude',...
                   'blocksToExclude','cutoff_zval','adaptdir'}));

% Defining input
addRequired(p,'trialInfo',@(x) validateattributes(x,{'numeric'},{'2d'}));
addRequired(p,'behavData',checkBehavData);
addRequired(p,'fileSpec',checkFileSpec);
% Parsing inputs
parse(p,trialInfo,behavData,fileSpec);

% Assigning input to variables
trialInfo = p.Results.trialInfo;
behavData = p.Results.behavData;
fileSpec = p.Results.fileSpec;

session = unique(trialInfo(:,1));
if numel(session) > 1
    error('buildtrialinfo:invalidInput',...
          'The input should contain data from one session only');
end
runs = unique(trialInfo(:,2));

temp = [behavData.run];
runIDs = [temp.id];

trialInfoOut = {};

for i = 1:numel(runs)
    actRun = runs(i);
    actBlocks = unique(trialInfo(trialInfo(:,2) == actRun,3));
    behavDataActRun = behavData(runIDs == actRun).data;
    temp = behavDataActRun(ismember(behavDataActRun.block,actBlocks),...
        {'session','run','block','resp','RT','aloc','vloc','blocktype','catch_trial'});
    temp.adaptdir = repmat({fileSpec.adaptdir},size(temp,1),1);
    temp.condition = trialInfo(trialInfo(:,2) == actRun & ismember(trialInfo(:,3),actBlocks),5);
    temp.hand = repmat(fileSpec.hand(fileSpec.runid == actRun),size(temp,1),1);
    temp.badtrials = trialInfo(trialInfo(:,2) == actRun & ismember(trialInfo(:,3),actBlocks),6);
    trialInfoOut = cat(1,trialInfoOut,mat2cell(dataset2struct(temp),ones(size(temp,1),1)));
end


end