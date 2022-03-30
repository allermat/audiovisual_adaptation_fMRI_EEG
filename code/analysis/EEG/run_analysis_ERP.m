saveDf = cd(fullfile(get_path('project'),'data'));
subjList = cellstr(ls);
subjList = regexp(subjList,'[0-9]{3}','match','once');
subjList = subjList(~ismember(subjList,{''}));
cd(saveDf);

for i = 1:numel(subjList)
    preproc_ERP_perside(subjList{i});
end
