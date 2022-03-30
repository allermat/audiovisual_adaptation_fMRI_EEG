function batch_glm_mvpa

% Load subject info
subjects = subject_info;
subjects = subjects([1:5]);

% Specify analysis
subanal = 'run_merged2';
if strfind(subanal, 'run')
    reg_unit = 'run';
else
    reg_unit = 'block';
end

for ss=1:numel(subjects)
    % SPM folder
    spmdir = fullfile(get_path('project'), get_folder(subjects(ss), 'r'), 'fMRI', '1st level', regexprep(subanal, '^_', '')); % s6wauf
    
    % Load onsets and regressor names
    fprintf('loading onsets...');
    if strfind(subanal, 'merged')
        [onset, names, blockstart] = load_onset(subjects(ss), reg_unit);
    else
        [onset, names] = load_onset(subjects(ss), reg_unit);
    end
    fprintf('done\n');
    
    %% ------- Initialize SPM and get functional files --------
    
    fprintf('SPM initialization...');
    spm_get_defaults('stats.maxmem', 2^33); % 8GB memory can be used at the same time during GLM estimation
    spm_get_defaults('stats.resmem', true); % use RAM (true) or hard disk (false)
    spm_jobman('initcfg');
    matlabbatch = struct('spm', {});
    fprintf('done\n');
    TR = 2.8;
    
    % Functional folders
    fprintf('get functional folders...');
    folders = {{'pretest'} {'posttest' 'ladapt'} {'posttest' 'radapt'}};
    funcfolder = cellfun(@(x) fullfile(get_path('project'), get_folder(subjects(ss), 'r', 'fMRI', 'processed data', x{:})), folders, 'un', 0);
    funcfolder = sort(cat(1, funcfolder{:})); % sorting is essential to match EPI files with onsets!!
    nsession = size(funcfolder, 1);
    if nsession ~= size(onset, 1)
        error('EPI folders do not match with onset data')
    end
    fprintf('done\n');
    
    % Functional and motion correction files
    fprintf('unzipping EPI files...');
    funcfile = get_3D_funcfile(funcfolder, 's6w*.nii'); % s6w s3
    fprintf('done\n');
    rpfile = cellfun(@(x) getfname(x, 'rp_f*.txt'), funcfolder);

    if strfind(subanal, 'merged')
        % Concatenate onsets and volumes of splitted runs
        issplitted = reshape(find(cellfun(@(x) ~isempty(strfind(x, 'adapttest')), funcfolder) & cellfun(@numel, blockstart) == 2), 2, []); % splitted runs
        nsplittedvols = cellfun(@numel, funcfile(issplitted(1,:)));
        for i=1:size(issplitted, 2)
            onset_restarted = cellfun(@(x) x+numel(funcfile{issplitted(1,i)})*TR, onset{issplitted(2,i)}, 'un', 0); % increase onset as if there were no scanner restart
%             onset{issplitted(1,i)} = cellfun(@(x,y) [x; y], onset{issplitted(1,i)}, onset_restarted, 'un', 0); % concatenate onsets
            switch reg_unit
                case 'block'
                    onset{issplitted(1,i)} = [onset{issplitted(1,i)} onset_restarted];
                    names{issplitted(1,i)} = [names{issplitted(1:2,i)}];
                case 'run'
                    onset{issplitted(1,i)} = cellfun(@(x,y) [x; y], onset{issplitted(1,i)}, onset_restarted, 'un', 0);
            end
            blockstart{issplitted(1,i)} = cat(1, blockstart{issplitted(1,i)}, blockstart{issplitted(2,i)}+numel(funcfile{issplitted(1,i)})*TR); % concatenate blockstart
            funcfile{issplitted(1,i)} = cat(1, funcfile{issplitted(1,i)}, funcfile{issplitted(2,i)}); % concatenate volumes
        end
        
        % Concatenate motion regressor files of splitted runs
        for i=1:size(issplitted, 2)
            rpdata = arrayfun(@(x) load(fullfile(funcfolder{x}, rpfile{x})), issplitted(:,i), 'un', 0);
            rpdata = cat(1, rpdata{:}); % concatenate data
            cnctfname = ['rp_cnct_' rpfile{issplitted(1,i)}(4:end-8) rpfile{issplitted(2,i)}(end-11:end-8) rpfile{issplitted(1,i)}(end-7:end)];
            fileID = fopen(fullfile(funcfolder{issplitted(1,i)}, cnctfname), 'w'); % open file for writing concatenated data
            for line=1:size(rpdata, 1)
                fprintf(fileID,' %e\t%e\t%e\t%e\t%e\t%e\n', rpdata(line,:));
            end
            fclose(fileID);
            rpfile{issplitted(1,i)} = cnctfname;
        end
        
        % Delete leftover data
        onset(issplitted(2,:)) = [];
        names(issplitted(2,:)) = [];
        blockstart(issplitted(2,:)) = [];
        funcfolder(issplitted(2,:)) = [];
        funcfile(issplitted(2,:)) = [];
        rpfile(issplitted(2,:)) = [];
        
        % Actialize sessions and volumes
        nsession = size(funcfolder, 1);
        nvolumes = cellfun(@numel, funcfile); % [0; cumsum(cellfun(@numel, funcfile))];
    end
    
    %% ----------- Specify model -------------
    
    id = 1;
    matlabbatch{id}.spm.stats.fmri_spec.dir{1} = spmdir;
    matlabbatch{id}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{id}.spm.stats.fmri_spec.timing.RT = TR;
    matlabbatch{id}.spm.stats.fmri_spec.timing.fmri_t = 38; % best to be the same as the number of slices
    matlabbatch{id}.spm.stats.fmri_spec.timing.fmri_t0 = 20; % should match the reference slice in slice time correction
    fprintf('specifying GLM\n');
    for s=1:nsession
        matlabbatch{id}.spm.stats.fmri_spec.sess(s).scans = funcfile{s}; % sid(s)
        for j=1:size(onset{s}, 2)
            matlabbatch{id}.spm.stats.fmri_spec.sess(s).cond(j).name = names{s}{j};
            matlabbatch{id}.spm.stats.fmri_spec.sess(s).cond(j).onset = onset{s}{j}; % sid(s)
            if strfind(names{s}{j}, 'vloc=')
                matlabbatch{id}.spm.stats.fmri_spec.sess(s).cond(j).duration = 10; % visual detection block
            else
                matlabbatch{id}.spm.stats.fmri_spec.sess(s).cond(j).duration = 0; % delta function
            end
            matlabbatch{id}.spm.stats.fmri_spec.sess(s).cond(j).tmod = 0;
            matlabbatch{id}.spm.stats.fmri_spec.sess(s).cond(j).pmod=struct('name',{},'param',{}, 'poly',{});
        end
        matlabbatch{id}.spm.stats.fmri_spec.sess(s).multi = {''};
        idsplittedrun = find(cellfun(@(x) ~isempty(strfind(x, 'cnct')), rpfile));
        if strfind(subanal, 'run_merged2') && ismember(s, idsplittedrun)
%             for b=2:numel(blockstart{s}) % add block specific regressors
                matlabbatch{id}.spm.stats.fmri_spec.sess(s).regress(1).name = 'Run1'; % b-1 sprintf('Block%d', b)
                regress = zeros(nvolumes(s), 1);
%                 if b == numel(blockstart{s})
%                     regress(floor(blockstart{s}(b)/TR):nvolumes(s)) = ones;
%                 else
%                     regress(floor(blockstart{s}(b)/TR):floor(blockstart{s}(b+1)/TR)-1) = ones;
%                 end
                regress(1:nsplittedvols(ismember(idsplittedrun, s))) = 1;
                matlabbatch{id}.spm.stats.fmri_spec.sess(s).regress(1).val = regress; % b-1
%             end
        else
            matlabbatch{id}.spm.stats.fmri_spec.sess(s).regress = struct('name',{},'val',{});
        end
        matlabbatch{id}.spm.stats.fmri_spec.sess(s).multi_reg = fullfile(funcfolder{s}, rpfile(s));
        matlabbatch{id}.spm.stats.fmri_spec.sess(s).hpf = 128;
    end
    matlabbatch{id}.spm.stats.fmri_spec.fact = struct('name',{},'levels',{});
    matlabbatch{id}.spm.stats.fmri_spec.bases.hrf.derivs = [1,0]; % time-derivative is included
    matlabbatch{id}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{id}.spm.stats.fmri_spec.global = 'None'; % no global normalization
    matlabbatch{id}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{id}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{id}.spm.stats.fmri_spec.cvi = 'AR(1)';
    spm_jobman('run', matlabbatch(id));
    fprintf('done\n');
   
    % ------- Model estimation ---------
    
    fprintf('estimating GLM\n');
    id = length(matlabbatch) + 1;
    spmfile = fullfile(spmdir, 'SPM.mat');
    matlabbatch{id}.spm.stats.fmri_est.spmmat{1} = spmfile;
    matlabbatch{id}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{id}.spm.stats.fmri_est.method.Classical = 1;
    spm_jobman('run', matlabbatch(id));
    fprintf('done\n');
    
%     delete_4D_funcfile(funcfolder, 's3*.nii'); % s6w s3
    
    % -------- Save data ---------
    
    str = datestr(now, 'yyyy_mm_dd_HH_MM');
    savedir = fullfile(get_path('project'), get_folder(subjects(ss), 'r'), 'fMRI', 'batch');
    mkdir(savedir);
    save(fullfile(savedir, sprintf('glm_%s_%s.mat', subanal, str)), 'matlabbatch');
end