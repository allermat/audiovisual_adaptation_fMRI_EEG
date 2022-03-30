function batch_inference

% Load subject info
subjects = subject_info;
subjects = subjects(1:5);
nsubjects = numel(subjects);

jobid = [1:2];

% ------- Calculate 1st level contrasts ---------

for ss=1:nsubjects
    % SPM file
    spmdir1{ss} = fullfile(get_path('project'), get_folder(subjects(ss), 'r'), ...
        'fMRI', '1st level', 'run_merged2'); % 's6wauf' 'run' run_merged run_s6wauf
    load(fullfile(spmdir1{ss}, 'SPM.mat'))
    matlabbatch{1}.spm.stats.con.spmmat{1} = fullfile(spmdir1{ss}, 'SPM.mat');
    
    % Session based substring filter
    adapt = {'radapt' 'ladapt'; 'radapt' 'ladapt'};
    session = subjects(ss).fMRI.pretest.session;
    for ses=1:numel(session)
        s = sum(ismember(subjects(ss).fMRI.pretest.match(1:ses), subjects(ss).fMRI.pretest.match{ses})); % session within adaptation 
        a = find(cellfun(@(x) x(1), upper(adapt(1,:))) == subjects(ss).fMRI.pretest.match{ses}); % adaptation side
        substr{s,a,1} = sprintf('.*session=%d,run=(%s).*', session(ses), strjoin(strsplit(num2str(subjects(ss).fMRI.pretest.runid{ses}), ' '), '|'));
    end
    for a=1:2 % adaptation side
        for ses=1:2
            substr{ses,a,2} = sprintf('.*session=%d,run=(%s).*', subjects(ss).fMRI.posttest.(adapt{ses,a}).session(ses), ...
                strjoin(strsplit(num2str(subjects(ss).fMRI.posttest.(adapt{ses,a}).runid{ses}), ' '), '|'));
        end
    end
%     substr = '';
 
    conid = 1;
    aloc = [-12 -5 -2 0 2 5 12];
    for i=aloc
        matlabbatch{1}.spm.stats.con.consess{conid}.tcon.name = sprintf('pre aloc%d resp > baseline', i);
        matlabbatch{1}.spm.stats.con.consess{conid}.tcon.convec = session_convec(SPM, ...
            cellfun(@(x) [x 'resp=1,aloc=' num2str(i) '\*bf\(1\)'], substr(:,:,1), 'un', 0));
        matlabbatch{1}.spm.stats.con.consess{conid}.tcon.sessrep = 'none';
        conid = conid + 1;
        matlabbatch{1}.spm.stats.con.consess{conid}.tcon.name = sprintf('pre aloc%d noresp > baseline', i);
        matlabbatch{1}.spm.stats.con.consess{conid}.tcon.convec = session_convec(SPM, ...
            cellfun(@(x) [x 'resp=0,aloc=' num2str(i) '\*bf\(1\)'], substr(:,:,1), 'un', 0));
        matlabbatch{1}.spm.stats.con.consess{conid}.tcon.sessrep = 'none';
        conid = conid + 1;
        matlabbatch{1}.spm.stats.con.consess{conid}.tcon.name = sprintf('rpost aloc%d resp > baseline', i);
        matlabbatch{1}.spm.stats.con.consess{conid}.tcon.convec = session_convec(SPM, ...
            cellfun(@(x) [x 'resp=1,aloc=' num2str(i) '\*bf\(1\)'], substr(:,1,2), 'un', 0));
        matlabbatch{1}.spm.stats.con.consess{conid}.tcon.sessrep = 'none';
        conid = conid + 1;
        matlabbatch{1}.spm.stats.con.consess{conid}.tcon.name = sprintf('rpost aloc%d noresp > baseline', i);
        matlabbatch{1}.spm.stats.con.consess{conid}.tcon.convec = session_convec(SPM, ...
            cellfun(@(x) [x 'resp=0,aloc=' num2str(i) '\*bf\(1\)'], substr(:,1,2), 'un', 0));
        matlabbatch{1}.spm.stats.con.consess{conid}.tcon.sessrep = 'none';
        conid = conid + 1;
        matlabbatch{1}.spm.stats.con.consess{conid}.tcon.name = sprintf('lpost aloc%d resp > baseline', i);
        matlabbatch{1}.spm.stats.con.consess{conid}.tcon.convec = session_convec(SPM, ...
            cellfun(@(x) [x 'resp=1,aloc=' num2str(i) '\*bf\(1\)'], substr(:,2,2), 'un', 0));
        matlabbatch{1}.spm.stats.con.consess{conid}.tcon.sessrep = 'none';
        conid = conid + 1;
        matlabbatch{1}.spm.stats.con.consess{conid}.tcon.name = sprintf('lpost aloc%d noresp > baseline', i);
        matlabbatch{1}.spm.stats.con.consess{conid}.tcon.convec = session_convec(SPM, ...
            cellfun(@(x) [x 'resp=0,aloc=' num2str(i) '\*bf\(1\)'], substr(:,2,2), 'un', 0));
        matlabbatch{1}.spm.stats.con.consess{conid}.tcon.sessrep = 'none';
        conid = conid + 1;
    end
    matlabbatch{1}.spm.stats.con.consess{conid}.tcon.name = 'aloc resp > baseline';
    matlabbatch{1}.spm.stats.con.consess{conid}.tcon.convec = session_convec(SPM, ...
        cellfun(@(x) [x 'resp=1,aloc=-?[0-9]+\*bf\(1\)'], substr(:,:,1:2), 'un', 0));
    matlabbatch{1}.spm.stats.con.consess{conid}.tcon.sessrep = 'none';
    conid = conid + 1;
    matlabbatch{1}.spm.stats.con.consess{conid}.tcon.name = 'aloc noresp > baseline';
    matlabbatch{1}.spm.stats.con.consess{conid}.tcon.convec = session_convec(SPM, ...
        cellfun(@(x) [x 'resp=0,aloc=-?[0-9]+\*bf\(1\)'], substr(:,:,1:2), 'un', 0));
    matlabbatch{1}.spm.stats.con.consess{conid}.tcon.sessrep = 'none';
    conid = conid + 1;
   
    matlabbatch{1}.spm.stats.con.delete = 1;
    
    % -------- Results --------
    
    con = matlabbatch{1}.spm.stats.con.consess;
    
    for i=1:numel(con)
        matlabbatch{2}.spm.stats.results.spmmat{1} = fullfile(spmdir1{ss}, 'SPM.mat');
        contype = fieldnames(con{i});
        matlabbatch{2}.spm.stats.results.conspec.titlestr = con{i}.(contype{:}).name;
        matlabbatch{2}.spm.stats.results.conspec.contrasts = i;
        matlabbatch{2}.spm.stats.results.conspec.threshdesc = 'none';
        matlabbatch{2}.spm.stats.results.conspec.thresh = 0.001;
        matlabbatch{2}.spm.stats.results.conspec.extent = 0;
        matlabbatch{2}.spm.stats.results.conspec.mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
        matlabbatch{2}.spm.stats.results.units = 1;
        matlabbatch{2}.spm.stats.results.print = false; % true
    end
    
    % -------- Run and save jobs ---------
    
    if any(ismember(jobid, [1 2]))
        % Run jobs
        arrayfun(@(x) spm_jobman('run', matlabbatch(x)), jobid(ismember(jobid, 1:2)));
        
        % Save jobs
        savedir = fullfile(get_path('project'), get_folder(subjects(ss), 'r'), 'fMRI', 'batch');
        matlabbatch = matlabbatch(jobid(ismember(jobid, 1:2)));
        save(fullfile(savedir, sprintf('con_%s.mat', datestr(now, 'yyyy_mm_dd_HH_MM'))), 'matlabbatch');
    end
    
    clear matlabbatch;
end


function convec = session_convec(SPM, expression)

convec = false(size(SPM.xX.name, 2), numel(expression));
for i=1:numel(expression)
    convec(:,i) = ~cellfun(@isempty, regexp(SPM.xX.name, expression{i}));
end
convec = double(any(convec, 2));

% Normalize contrast
convec = convec / sum(convec);