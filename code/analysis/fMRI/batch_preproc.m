function batch_preproc(varargin)

oldpreproc = 0;

% Load subject info
subjects = subject_info;
subjects = subjects([1:5]);

varargin = {'realign' 'temporal' 'segment' 'smooth'};

for ss=1%:length(subjects)
    % Set up anatomical
    anatfolder = fullfile(get_path('project'), get_folder(subjects(ss), 'r', 'anat'));
    anatfile = getfname(anatfolder, 's*.nii');

    % Set up functional folders
    folders = {{'pretest'} {'posttest' 'radapt'} {'posttest' 'ladapt'}};  % 
    funcfolder = cellfun(@(x) fullfile(get_path('project'), get_folder(subjects(ss), 'r', 'fMRI', 'processed data', x{:})), folders, 'un', 0);
    funcfolder = sort(cat(1, funcfolder{:}));
    nsessions = size(funcfolder, 1);
%     meanfuncfolder = fullfile(get_path('project'), get_folder(subjects(ss), 'r'), 'fMRI', 'processed data');
    
    % SPM and matlabbatch initialization
    spm('defaults', 'fmri');
    spm_jobman('initcfg');
%     spm_figure('GetWin','Graphics'); % make sure the post script files are saved
    matlabbatch = struct('spm', {});
    
    % --------- Motion correction and unwarping ---------
    if ismember('realign', varargin)
        id = length(matlabbatch) + 1;
        funcfile = get_3D_funcfile(funcfolder, 'f*.nii.gz', subjects(ss));
        for s=1:nsessions
%             if oldpreproc && s == 1
%                 matlabbatch{id}.spm.spatial.realignunwarp.data(s).scans = funcfile{s}(1); % reference image
%             else
            matlabbatch{id}.spm.spatial.realignunwarp.data(s).scans = funcfile{s};
%             end
            matlabbatch{id}.spm.spatial.realignunwarp.data(s).pmscan = '';
        end
        matlabbatch{id}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
        matlabbatch{id}.spm.spatial.realignunwarp.eoptions.sep = 4;
        matlabbatch{id}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
        matlabbatch{id}.spm.spatial.realignunwarp.eoptions.rtm = 0; % register to first image!!
        matlabbatch{id}.spm.spatial.realignunwarp.eoptions.einterp = 2;
        matlabbatch{id}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
        matlabbatch{id}.spm.spatial.realignunwarp.eoptions.weight = '';
        matlabbatch{id}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
        matlabbatch{id}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
        matlabbatch{id}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
        matlabbatch{id}.spm.spatial.realignunwarp.uweoptions.jm = 0;
        matlabbatch{id}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
        matlabbatch{id}.spm.spatial.realignunwarp.uweoptions.sot = [];
        matlabbatch{id}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
        matlabbatch{id}.spm.spatial.realignunwarp.uweoptions.rem = 1;
        matlabbatch{id}.spm.spatial.realignunwarp.uweoptions.noi = 5;
        matlabbatch{id}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
        matlabbatch{id}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1]; % all images + mean image
        matlabbatch{id}.spm.spatial.realignunwarp.uwroptions.rinterp = 4; % 4th degree B-spline interpolation
        matlabbatch{id}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
        matlabbatch{id}.spm.spatial.realignunwarp.uwroptions.mask = 1;
        matlabbatch{id}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';
        spm_jobman('run', matlabbatch(id));
        gzip_4D_funcfile(funcfolder, 'uf*.nii');
        delete_4D_funcfile(funcfolder, 'f*.nii');
        if oldpreproc
           funcfolder(1) = [];
           nsessions = nsessions - 1;
        end
%         meanfunc = fullfile(funcfolder{1}, getfname(funcfolder{1}, 'meanuf*.nii'));
%         movefile(meanfunc{1}, meanfuncfolder);
    end
    
    % ------------ Slice time correction ------------
    if ismember('temporal', varargin)
        id = length(matlabbatch) + 1;
        funcfile = get_3D_funcfile(funcfolder, 'uf*.nii.gz');
        for s=1:nsessions
            matlabbatch{id}.spm.temporal.st.scans{s} = funcfile{s};
        end
        matlabbatch{id}.spm.temporal.st.tr = 2.8;
        matlabbatch{id}.spm.temporal.st.nslices = 38;
        matlabbatch{id}.spm.temporal.st.ta = 2.8-2.8/matlabbatch{id}.spm.temporal.st.nslices;
        matlabbatch{id}.spm.temporal.st.so = 1:matlabbatch{id}.spm.temporal.st.nslices; % ascending acquisition order
        matlabbatch{id}.spm.temporal.st.refslice = 20; % reference slice is the middle one
        matlabbatch{id}.spm.temporal.st.prefix = 'a';
        spm_jobman('run', matlabbatch(id));
        gzip_4D_funcfile(funcfolder, 'auf*.nii');
        delete_4D_funcfile(funcfolder, 'uf*.nii');
    end
    
    % ------------ Segmentation and skull stripping ------------
    if ismember('segment', varargin)
        id = length(matlabbatch) + 1;
        matlabbatch{id}.spm.spatial.preproc.channel.vols = fullfile(anatfolder, anatfile);
        matlabbatch{id}.spm.spatial.preproc.channel.biasreg = 0.001;
        matlabbatch{id}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{id}.spm.spatial.preproc.channel.write = [0 1]; % write out bias-corrected image
        ngaus = [1 1 2 3 4 2];
        native = {[1 0] [1 0] [0 0] [0 0] [0 0] [0 0]}; % write out grey matter (1) and white matter (2)
        warped = {[0 0] [0 0] [0 0] [0 0] [0 0] [0 0]};
        for t=1:6
            matlabbatch{id}.spm.spatial.preproc.tissue(t).tpm = fullfile(spm('dir'), 'tpm', {sprintf('TPM.nii,%d', t)});
            matlabbatch{id}.spm.spatial.preproc.tissue(t).ngaus = ngaus(t);
            matlabbatch{id}.spm.spatial.preproc.tissue(t).native = native{t};
            matlabbatch{id}.spm.spatial.preproc.tissue(t).warped = warped{t};
        end
        matlabbatch{id}.spm.spatial.preproc.warp.mrf = 1;
        matlabbatch{id}.spm.spatial.preproc.warp.cleanup = 1;
        matlabbatch{id}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{id}.spm.spatial.preproc.warp.affreg = 'mni'; % 'mni' for European brains, could be either 'eastern' or 'subj'
        matlabbatch{id}.spm.spatial.preproc.warp.fwhm = 0;
        matlabbatch{id}.spm.spatial.preproc.warp.samp = 3;
        matlabbatch{id}.spm.spatial.preproc.warp.write = [1 1]; % write deformation field in both direction
        spm_jobman('run', matlabbatch(id));
        id = length(matlabbatch) + 1;
        bcfile = getfname(anatfolder, 'ms*.nii'); % bias corrected image
        gmfile = getfname(anatfolder, 'c1s*.nii'); % grey matter image
        wmfile = getfname(anatfolder, 'c2s*.nii'); % white matter image
        matlabbatch{id}.spm.util.imcalc.input = fullfile(anatfolder, [bcfile; gmfile; wmfile]);
        matlabbatch{id}.spm.util.imcalc.output = ['c12' bcfile{1}]; % output file with pattern replacement
        matlabbatch{id}.spm.util.imcalc.outdir{1} = anatfolder;
        matlabbatch{id}.spm.util.imcalc.expression = 'i1 .* ((i2 + i3) > 0)'; % thresholding
        matlabbatch{id}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{id}.spm.util.imcalc.options.mask = 0;
        matlabbatch{id}.spm.util.imcalc.options.interp = 1;
        matlabbatch{id}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('run', matlabbatch(id)); 
    end
    
    % ---------- Smoothing ----------
    if ismember('smooth', varargin)
        id = length(matlabbatch) + 1;
        try delete_4D_funcfile(funcfolder, 'auf*.mat'), catch, end;
%         continue
        funcfile = get_3D_funcfile(funcfolder, 'auf*.nii.gz'); % wrauf auf
%         meanfunc = fullfile(meanfuncfolder, getfname(funcfolder{1}, 'meanuf*.nii'));
        matlabbatch{id}.spm.spatial.smooth.data = cat(1, funcfile{:});
        matlabbatch{id}.spm.spatial.smooth.fwhm = [3 3 3]; % 6 6 6 3 3 3
        matlabbatch{id}.spm.spatial.smooth.dtype = 0;
        matlabbatch{id}.spm.spatial.smooth.im = 0;
        matlabbatch{id}.spm.spatial.smooth.prefix = 's3'; % s3 s6
        spm_jobman('run', matlabbatch(id));
        gzip_4D_funcfile(funcfolder, 's3auf*.nii'); % s6wrauf s3auf
        delete_4D_funcfile(funcfolder, 'auf*.nii'); % wrauf auf
    end
end