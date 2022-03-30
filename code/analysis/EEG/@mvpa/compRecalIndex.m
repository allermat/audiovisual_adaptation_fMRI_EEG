function [RI,RIbin] = compRecalIndex(predLabelInfo,predLabels)
% Method to compute recalibration index

idx_post_L = predLabelInfo.blocktype == 'post-test' & predLabelInfo.adaptdir == 'L';
idx_post_R = predLabelInfo.blocktype == 'post-test' & predLabelInfo.adaptdir == 'R';
idx_aloc = predLabelInfo.aloc;
aloc_levels = unique(predLabelInfo.aloc(~isnan(predLabelInfo.aloc)));
predLabelsBin = predLabels >= 0;

if all(~idx_post_L) || all(~idx_post_R)
    warning('mvpa:compRecalIndex:requestedDataNotPresent',...
        ['Left or right recalibration data is missing ',...
        'recalibration index can''t be computed, returning.']);
    [RI,RIbin] = deal([]);
    return;
end

mean_perloc_post_L = arrayfun(@(x) shiftdim(mean(predLabels(idx_post_L & idx_aloc == x,:,:)),1),...
                              aloc_levels,'UniformOutput',false);
mean_perloc_post_R = arrayfun(@(x) shiftdim(mean(predLabels(idx_post_R & idx_aloc == x,:,:)),1),...
                              aloc_levels,'UniformOutput',false);
% mean_post_L = shiftdim(mean(predLabels(idx_post_L,:,:)),1);
% mean_post_R = shiftdim(mean(predLabels(idx_post_R,:,:)),1);
% var_post_L = shiftdim(var(predLabels(idx_post_L,:,:)),0,1);
% var_post_R = shiftdim(var(predLabels(idx_post_R,:,:)),0,1);
mean_post_L = mean(cat(3,mean_perloc_post_L{:}),3);
mean_post_R = mean(cat(3,mean_perloc_post_R{:}),3);
var_post_L = var(cat(3,mean_perloc_post_L{:}),0,3);
var_post_R = var(cat(3,mean_perloc_post_R{:}),0,3);

% Computing the percentage of trials predicted to right for each
% location and adaptation direction
fun = @(x) shiftdim(sum(predLabelsBin(idx_post_L & idx_aloc == x,:,:)),1)/...
      sum(idx_post_L & idx_aloc == x)*100;
bin_perloc_post_L = arrayfun(fun,aloc_levels,'UniformOutput',false);
fun = @(x) shiftdim(sum(predLabelsBin(idx_post_R & idx_aloc == x,:,:)),1)/...
      sum(idx_post_R & idx_aloc == x)*100;
bin_perloc_post_R = arrayfun(fun,aloc_levels,'UniformOutput',false);

pct_post_L = mean(cat(3,bin_perloc_post_L{:}),3);
pct_post_R = mean(cat(3,bin_perloc_post_R{:}),3);
var_pct_post_L = var(cat(3,bin_perloc_post_L{:}),0,3);
var_pct_post_R = var(cat(3,bin_perloc_post_R{:}),0,3);


RI = mean_post_R - mean_post_L;
% RI_err = sqrt(var_post_R/sum(idx_post_R) + var_post_L/sum(idx_post_L));
RIbin = pct_post_R - pct_post_L;

end