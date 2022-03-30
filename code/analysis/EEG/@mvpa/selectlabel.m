function labelOut = selectlabel(inputLabel,inputConds)
% Selects the label for the given input label and condition
% 
% INPUT:
%   inputLabel = input label
%   inputCond = input condition
% 
% OUTPUT:
%   labelOut = the output label (for training or generalization)
% 

% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

%% Parsing input
p = inputParser;

validLabels = {'hand','aloc','vloc'};
validConds = ['adapt|pre|post|adL|adR|hL|hR|catch|nCatch|',...
              '(ses[0-9]+,?)|(aloc[0-9n]+,?)'];

addRequired(p,'inputLabel',@(x)any(validatestring(x,validLabels)));
addRequired(p,'inputConds',...
            @(x)all(~cellfun(@isempty,regexp(regexp(x,'-','split'),validConds,'once'))));

% Parsing inputs.
parse(p,inputLabel,inputConds);

% Assigning inputs to variables
inputLabel = p.Results.inputLabel;
inputConds = p.Results.inputConds;

%% Main
% Splitting the input conditions to parts
conds = regexp(inputConds,'-','split');

if strcmp(inputLabel,'aloc')
    labelOut = 'aloc';
elseif strcmp(inputLabel,'vloc')
    labelOut = 'vloc';
elseif strcmp(inputLabel,'resp')
    labelOut = 'resp';
elseif strcmp(inputLabel,'hand')
    labelOut = 'hand';
else
    error('Unidentified label');
end

end