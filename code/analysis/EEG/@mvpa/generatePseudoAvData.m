function [info,infoPseu,featPseu] = generatePseudoAvData(info,feat,condDef)

% Selecting A and V only trials
isAxorV = isnan(info.locV) | isnan(info.locA);
info = info(isAxorV,:);
feat = feat(:,:,isAxorV);

locAlevels = unique(info.locA(~isnan(info.locA)));
locVlevels = unique(info.locV(~isnan(info.locV)));

% Pre-allocating arrays
[featPseuTemp,infoPseuTemp] = deal(cell(length(locVlevels),length(locAlevels)));

for i = 1:length(locVlevels)
    for j = 1:length(locAlevels)
        
        trialIdxV = find(info.locV == locVlevels(i));
        trialIdxA = find(info.locA == locAlevels(j));
        
        if length(trialIdxA) < length(trialIdxV)
            trialIdxV = randsample(trialIdxV,length(trialIdxA));
        else
            trialIdxA = randsample(trialIdxA,length(trialIdxV));
        end
        
        featPseuTemp{i,j} = sum(cat(4,feat(:,:,trialIdxA),feat(:,:,trialIdxV)),4);
        tempInfo = info(trialIdxV,{'condition','hand','task','locA','locV','relV','resp'});
        tempInfo.locA = info.locA(trialIdxA);
        infoPseuTemp{i,j} = tempInfo;
        
    end
end

infoPseu = cat(1,infoPseuTemp{:});
featPseu = cat(3,featPseuTemp{:});

infoPseu.task = NaN(size(infoPseu.task));
infoPseu.hand = NaN(size(infoPseu.hand));
infoPseu.resp = NaN(size(infoPseu.resp));
for i = 1:size(infoPseu,1)
    % Assigning condition labels to the pseudoAV trials. I use the
    % condition labels corresponding to the locA, locV and relV and keep
    % the task as auditory localization as in the generalization this info 
    % will not be specified. This is a workaround, might require a proper
    % solution. 
    infoPseu.condition(i) = condDef.condition(...
        condDef.task == 1 & ...
        condDef.locationAuditory == infoPseu.locA(i) & ...
        condDef.locationVisual == infoPseu.locV(i) & ...
        condDef.reliabilityVisual == infoPseu.relV(i));
end
infoPseu.example = (1:size(infoPseu,1))';
infoPseu = infoPseu(:,[end,1:end-1]);

end

