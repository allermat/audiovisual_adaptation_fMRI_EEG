function [info,infoCorr,featCorr] = generateCorrectedAvData(info,feat,condDef,corrMode)
% Method for generating corrected AV data

% Selecting A and V only trials for averaging
condsAxorV = condDef.condition(isnan(condDef.locationAuditory) | ...
                               isnan(condDef.locationVisual));
temp = info(ismember(info.condition,condsAxorV),{'condition', ...
                    'hand','task','locA','locV','relV'});
[~,idx] = unique(temp(:,{'condition','hand'}),'rows');
infoAvg = temp(idx,:);
featAvg = NaN(size(feat,1),size(feat,2),size(infoAvg,1));

for i = 1:size(infoAvg,1)
    actTrials = info.condition == infoAvg.condition(i) & ...
        info.hand == infoAvg.hand(i);
    featAvg(:,:,i) = mean(feat(:,:,actTrials),3);        
end

% Only AV trials will be included in the final dataset
isAVtrial = ismember(info.condition,condDef.condition(~isnan(condDef.locationAuditory) & ...
                               ~isnan(condDef.locationVisual)));
infoCorr = info(isAVtrial,{'session','iTrialSession','condition', ...
                    'hand','task','locA','locV','relV','resp'});
featCorr = feat(:,:,isAVtrial);
for i = 1:size(infoCorr,1)
    actAudAvgIdx = infoAvg.locA == infoCorr.locA(i) & ...
        infoAvg.hand == infoCorr.hand(i);
    actVisAvgIdx = infoAvg.locV == infoCorr.locV(i) & ...
        infoAvg.hand == infoCorr.hand(i) & ...
        infoAvg.relV == infoCorr.relV(i);
    % Correction mode 1: AV - (A + V)
    if corrMode == 1
        featCorr(:,:,i) = featCorr(:,:,i)-featAvg(:,:,actAudAvgIdx);
    end
    % Correction mode 2: AV - V
    featCorr(:,:,i) = featCorr(:,:,i)-featAvg(:,:,actVisAvgIdx);
end

end

