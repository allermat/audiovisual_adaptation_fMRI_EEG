function ROI_imcalc

% Load subject info
subjects = subject_info;
subjects = subjects(1:5);

% SPM and batch initialization
spm('defaults', 'fmri');
spm_jobman('initcfg');
matlabbatch = struct('spm', {});

for ss=1:length(subjects)
    filestart = 'rSROI_l'; % define hemisphere
    basefolder = fullfile(get_path('project'), get_folder(subjects(ss), 'r'), 'fMRI', 'ROI');
    id = length(matlabbatch) + 1;
    fname = {'h.IPS0_maxprob'; 'h.IPS1_maxprob';  'h.IPS2_maxprob'; ...
        'h.IPS3_maxprob'; 'h.IPS4_maxprob'; 'h.IPS5_maxprob'; 'h.SPL1_maxprob'};
    matlabbatch{id}.spm.util.imcalc.input = cellfun(@(x) fullfile(basefolder, ...
        'ProbAtlas_v4', 'freesurfer10', 'coregistered', [filestart x '.nii']), fname, 'un', 0);
    matlabbatch{id}.spm.util.imcalc.output = [filestart 'h.IPS0-5_maxprob.nii'];
    matlabbatch{id}.spm.util.imcalc.outdir{1} = fullfile(basefolder, ...
        'ProbAtlas_v4', 'freesurfer10', 'coregistered');
    matlabbatch{id}.spm.util.imcalc.expression = 'i1 | i2 | i3 | i4 | i5 | i6 | i7';
    matlabbatch{id}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{id}.spm.util.imcalc.options.mask = 0;
    matlabbatch{id}.spm.util.imcalc.options.interp = 1;
    matlabbatch{id}.spm.util.imcalc.options.dtype = 4;
    spm_jobman('run', matlabbatch(id));
end









