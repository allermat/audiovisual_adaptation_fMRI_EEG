function NMF = fitNeurometricFunctions(predLabels,infoGenExamples,varargin)
% Method to fit neurometric functions across training and
% generalization samples.

% Parsing input
validGenTimes = {'tr','tr_x_tr'};
p = inputParser;
addRequired(p,'predLabels',@(x) validateattributes(x,{'numeric'}, ...
                                                  {'nonempty'}));
addRequired(p,'infoGenExamples',@(x) validateattributes(x,{'table'}, ...
                                                  {'nrows',size(predLabels,1)}));
addOptional(p,'genTimeStr','tr',@(x) any(validatestring(x,validGenTimes)));
parse(p,predLabels,infoGenExamples,varargin{:});
predLabels = p.Results.predLabels;
infoGenExamples = p.Results.infoGenExamples;
genTimeStr = p.Results.genTimeStr;

% Binarizing predicted labels
predRight = predLabels > 0;

% feeding binarized predicted labels into the info table to utilize
% the functionalities of varfun
if strcmp(genTimeStr,'tr')
    s = {ones(size(predRight,1),1),size(predRight,2)};
else
    s = {ones(size(predRight,1),1),size(predRight,2), ...
         size(predRight,3)};
end
infoGenExamples.predRight = mat2cell(predRight,s{:});

% Generating groupings for the conditions of the neurometric curves
infoGenExamples.nmGrouping = cell(size(infoGenExamples,1),1);
infoGenExamples.nmGrouping(infoGenExamples.blocktype == 'pre-test') = ...
    repmat({'pre'},sum(infoGenExamples.blocktype == 'pre-test'),1);
infoGenExamples.nmGrouping(infoGenExamples.blocktype == 'post-test'& ...
                           infoGenExamples.adaptdir == 'R') = ...
    repmat({'post_R'},sum(infoGenExamples.blocktype == 'post-test'& ...
                     infoGenExamples.adaptdir == 'R'),1);
infoGenExamples.nmGrouping(infoGenExamples.blocktype == 'post-test'& ...
                           infoGenExamples.adaptdir == 'L') = ...
    repmat({'post_L'},sum(infoGenExamples.blocktype == 'post-test'& ...
                     infoGenExamples.adaptdir == 'L'),1);
infoGenExamples.nmGrouping = categorical(infoGenExamples.nmGrouping);

num = varfun(@(x) sum(cat(1,x{:})),infoGenExamples,'InputVariables',{'predRight'},...
           'GroupingVariables',{'nmGrouping','aloc'},'OutputFormat','cell');
outOfNum = varfun(@(x) size(x,1),infoGenExamples,'InputVariables',{'predRight'},...
           'GroupingVariables',{'nmGrouping','aloc'});
stimLevels = unique(infoGenExamples.aloc);
nmGroupingLevels = categories(outOfNum.nmGrouping);

% We fit the cumulative normal function to the data
PF = @PAL_CumulativeNormal;
% Fixed parameters of the fit
% guess rate
opt.gamma = 0.2;
% lapse rate
opt.lambda = 0.2;

% Starting the timer and printing details
cStart = clock;
ds = datestr(now);
printFnTitle(80,'fitNeurometricFunctions',ds)
fprintf('Estimating... \n');

NMF = struct();
NMF.PF = PF;
NMF.stimLevels = stimLevels;
for i = 1:numel(nmGroupingLevels)
    actNum = num(outOfNum.nmGrouping == nmGroupingLevels{i});
    actNum = cat(1,actNum{:});
    actOutOfNum = outOfNum.GroupCount(outOfNum.nmGrouping == nmGroupingLevels{i});
    [pv,pctr] = mvpa.fitNMFacrossTime(PF,stimLevels,actNum,actOutOfNum,opt);
    NMF.([nmGroupingLevels{i},'_pv']) = pv;
    NMF.([nmGroupingLevels{i},'_pctr']) = pctr;
end

% Finishing timer and printing elapsed time
fprintf('Estimation elapsed time (days hours:minutes:seconds) %s \n\n',...
    datestr(etime(clock,cStart)/86400,'dd HH:MM:SS'));

end