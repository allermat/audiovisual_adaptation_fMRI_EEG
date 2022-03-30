function isExample = selectexamples(inputConds,condDef,dataInfo)
% Selects examples matching the given input condition
% 
% INPUT:
%   cond: hyphen separated string of the input conditions
%   condDef: table of stimulus condition definitions
%   dataInfo: table of information about the examples
% 
% OUTPUT:
%   isExample: N x 1 logical vector indicating which examples correspond to
%       the given input conditions. N = number of examples in dataInfo.
%

% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

%% Parsing input
p = inputParser;

validConds = ['adapt|all|pre|post|adL|adR|hL|hR|catch|nCatch|',...
              '(ses[0-9]+,?)|(aloc[0-9n]+,?)'];

addRequired(p,'inputConds',...
            @(x)all(~cellfun(@isempty,regexp(regexp(x,'-','split'),validConds,'once'))));
addRequired(p,'condDef',@(x)validateattributes(x,{'table'},{'nonempty'}));
addRequired(p,'dataInfo',@(x)validateattributes(x,{'table'},{'nonempty'}));

% Parsing inputs
parse(p,inputConds,condDef,dataInfo);

% Assigning inputs to variables
inputConds = p.Results.inputConds;
condDef = p.Results.condDef;
dataInfo = p.Results.dataInfo;

%% Main
% Splitting the input conditions to parts
condsCell = regexp(inputConds,'-','split');
% Making sure, that only one level of a condition category is specified at 
% a time
moreThanOnePerCat = sum(ismember({'all','adapt','pre','post'},condsCell)) > 1 || ...
                    sum(ismember({'all','adL','adR'},condsCell)) > 1 || ...
                    sum(ismember({'all','hL','hR'},condsCell)) > 1 || ...
                    sum(ismember({'all','catch','nCatch'},condsCell)) > 1;
if moreThanOnePerCat
    error('mvpa:selectexamples:ambiguousCondition',...
        'Only one level of a condition category can be specified at a time');
end

% Automatically excluding catch trials
isExample = true(size(dataInfo,1),1);

for i = 1:size(condsCell,2)
    if ismember(condsCell(i),'all')
        % In this case do nothing, all examples are included
    elseif ismember(condsCell(i),'pre')
        isExample = isExample & dataInfo.blocktype == 'pre-test';
    elseif ismember(condsCell(i),'post')
        isExample = isExample & dataInfo.blocktype == 'post-test';
    elseif ismember(condsCell(i),'adapt')
        isExample = isExample & ...
            (dataInfo.blocktype == 'l-adaptation' | ...
             dataInfo.blocktype == 'r-adaptation');
    elseif ismember(condsCell(i),'adL')
        isExample = isExample & dataInfo.adaptdir == 'L';
    elseif ismember(condsCell(i),'adR')
        isExample = isExample & dataInfo.adaptdir == 'R';
    elseif ismember(condsCell(i),'hL')
        isExample = isExample & dataInfo.hand == 'L';
    elseif ismember(condsCell(i),'hR')
        isExample = isExample & dataInfo.hand == 'R';
    elseif ismember(condsCell(i),'catch')
        isExample = isExample & dataInfo.catch_trial;
    elseif ismember(condsCell(i),'nCatch')
        isExample = isExample & ~dataInfo.catch_trial;
    elseif strfind(condsCell(i),'ses')
        s = regexp(regexp(condsCell(i),'ses(\d+,?)+','tokens','once'), ...
                   ',','split');
        s = cellfun(@str2num,s{:});
        isExample = isExample & ismember(dataInfo.session,s);
    elseif strfind(condsCell(i),'aloc')
        loc = regexp(regexp(condsCell(i),'aloc([0-9n]+,?)+','tokens','once'), ...
                   ',','split');
        loc = cellfun(@str2num,strrep(loc{:},'n','-'));
        isExample = isExample & ismember(dataInfo.aloc,loc);
    else
        error('Unrecognized condition');
    end
    
end

end