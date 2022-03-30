clearvars;
modality = {'PP','fMRI','EEG'};
for i = 1:numel(modality)
    out = load([modality{i},'_behav_all.mat']);
    out = out.out;
    temp = cat(1,out.pre);
    temp.false_alarm;
    fa_pre = cat(1,temp.false_alarm);
    temp = cat(1,out.post);
    fa_post = cat(1,temp.false_alarm);
    fa_post = cat(1,cat(1,fa_post.ladapt),cat(1,fa_post.radapt));
    temp = cat(1,out.adapt);
    fa_adapt = cat(1,temp.false_alarm);
    fa_adapt = cat(1,cat(1,fa_adapt.ladapt),cat(1,fa_adapt.radapt));
    temp_fa = cat(1,fa_pre,fa_post,fa_adapt);
    temp_mod = repmat(lower(modality(i)),numel(temp_fa),1);
    temp_phase = repmat({'pre','postAV','postVA','adaptAV','adaptVA'},size(fa_pre,1),1);
    temp_phase = cat(1,temp_phase(:));
    fa_cell{i} = table(temp_mod,temp_phase,temp_fa,'VariableNames',...
                       {'modality','phase','FA_rate'});
end
fa_table = cat(1,fa_cell{:});
group_fa_table = varfun(@mean,fa_table,'InputVariables',{'FA_rate'},'GroupingVariables',...
                        {'modality','phase'});
group_fa_table = join(group_fa_table,...
                      varfun(@std,fa_table,'InputVariables',{'FA_rate'},'GroupingVariables',...
                            {'modality','phase'}));
group_fa_table.sem_FA_rate = group_fa_table.std_FA_rate./sqrt(group_fa_table.GroupCount);