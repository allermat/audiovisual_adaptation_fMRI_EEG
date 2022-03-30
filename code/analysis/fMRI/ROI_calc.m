function ROI_calc
%   ROI_calc does the following operations on ROIs using SPM: coregistration 
%   and merging. The operations are attached to jobs using particular 
%   atlases and ROIs.
%
% Although, coregistration might be performed multiple times, but using jobs
% allow independent operations on specific atlases and ROIs.
%
% Comments:
% when merging ROI columns are defined as follows:
%   - column 1-2: files to be merged
%   - column 3: patterns for string replacement in file1-output format

%% Initialization

% Define jobs
job = struct('atlas', {}, 'operation', {}, 'ROI', {});
job(1) = struct('atlas', 'FS_Destrieux', 'operation', {{'coreg'}}, ...
    'ROI', {{'*IPL.nii*'; '*hAud*'; '*G_temp_sup-G_T_transv*'}});
job(2) = struct('atlas', 'FS_Destrieux', 'operation', {{'merge'}}, ...
    'ROI', {{'rsROI_lh.IPL.nii' 'rsROI_rh.IPL.nii' 'lh.-'; ...
    'rsROI_lh.hAud*' 'rsROI_rh.hAud*' 'lh.-'; ...
    'rsROI_lh.G_temp_sup-G_T_transv*' 'rsROI_rh.G_temp_sup-G_T_transv*' 'lh.-'}});
job(3) = struct('atlas', fullfile('ProbAtlas_v4', 'freesurfer10'), 'operation', {{'coreg'}}, ...
    'ROI', {{'*maxprob*'}});
job(4) = struct('atlas', fullfile('ProbAtlas_v4', 'freesurfer10'), 'operation', {{'merge'}}, ...
    'ROI', {{'rSROI_lh*maxprob.nii' 'rSROI_rh*maxprob.nii' 'lh.-'; ...
    'rSROI_IPS0_m*' 'rSROI_IPS1_m*' '0-01'; ...
    'rSROI_IPS01_m*' 'rSROI_IPS2_m*' '01-02'; ...
    'rSROI_IPS2_m*' 'rSROI_IPS3_m*' '2-23'; ...
    'rSROI_IPS23_m*' 'rSROI_IPS4_m*' '23-24'; ...
    'rSROI_IPS24_m*' 'rSROI_IPS5_m*' '24-25'; ...
    'rSROI_IPS25_m*' 'rSROI_SPL1_m*' '25-25+'; ...
    'rSROI_IPS3_m*' 'rSROI_IPS4_m*' '3-34'; ...
    'rSROI_IPS34_m*' 'rSROI_IPS5_m*' '34-35'; ...
    'rSROI_IPS35_m*' 'rSROI_SPL1_m*' '35-35+'}});
job(5) = struct('atlas', fullfile('ProbAtlas_v4', 'freesurfer10'), 'operation', {{'merge'}}, ...
    'ROI', {{'rSROI_lh*maxprob.nii' 'rSROI_rh*maxprob.nii' 'lh.-';}});

% Load subject info
subjects = subject_info;
subjects = subjects([1:5]); % 1:4 6:8 10

% SPM and batch initialization
spm('defaults', 'fmri');
spm_jobman('initcfg');
matlabbatch = struct('spm', {});

for ss=1:length(subjects)
    subject = subjects(ss);
    for j=1:length(job)
        
        % Number of ROI patterns in the current job
        nROIs = size(job(j).ROI, 1);
        
        %% Coregistration
        if ismember('coreg', job(j).operation)
            % ROI folder
            ROIfolder = fullfile(get_path('project'), get_folder(subject, 'w'), ...
                'fMRI', 'ROI', job(j).atlas, 'coregistered');
            mkdir(ROIfolder);
            
            % Setup anatomical and EPI file to be coregistered to
            anatfolder = fullfile(get_path('project'), get_folder(subject, 'r', ...
                'anat'));
            anatfile = getfname(anatfolder, 'c12ms*.nii');
            infolder = fullfile(get_path('project'), get_folder(subject, 'r'), ...
                'fMRI', 'ROI');
            copyfile(fullfile(anatfolder, anatfile{:}), infolder);
            funcfolder = fullfile(get_path('project'), get_folder(subject, 'r', ...
                'fMRI', 'processed data', 'pretest'));
            meanfuncfile = getfname(funcfolder{1}, 'mean*');
            funcfile = get_3D_funcfile(funcfolder(1), 'uf*.nii.gz'); % this is in the same space as the mean functional...
            firstfuncfile = funcfile{1}(1);
            
            % Setup ROI files
            outfolder = fullfile(get_path('project'), get_folder(subject, 'r'), ...
                'fMRI', 'ROI', job(j).atlas);
            allfile = cell(nROIs, 1);
            for r=1:nROIs
                copyfile(fullfile(outfolder, strrep(job(j).ROI{r}, '_MNI', '')), ROIfolder);
                allfile{r} = getfname(ROIfolder, strrep(job(j).ROI{r}, '_MNI', ''));
                cellfun(@(x) movefile(fullfile(ROIfolder, x), fullfile(ROIfolder, ...
                    ['SROI_' x])), allfile{r}(cellfun(@isempty, strfind(allfile{r}, 'SROI'))), 'un', 1);
                allfile{r} = getfname(ROIfolder, ['SROI_' strrep(job(j).ROI{r}, '_MNI', '')]);
            end
            allfile = cat(1, allfile{:});
            
            % Coregister and reslice images
            id = length(matlabbatch) + 1;
            matlabbatch{id}.spm.spatial.coreg.estwrite.source(1) = ...
                fullfile(infolder, anatfile); % anatomy is coregistered to functional
            matlabbatch{id}.spm.spatial.coreg.estwrite.ref(1) = ...
                fullfile(funcfolder{1}, meanfuncfile); 
            matlabbatch{id}.spm.spatial.coreg.estwrite.other = ...
                fullfile(ROIfolder, allfile);
            matlabbatch{id}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi'; % normalized mutual information
            matlabbatch{id}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
            matlabbatch{id}.spm.spatial.coreg.estwrite.eoptions.tol = ...
                [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
            matlabbatch{id}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
            matlabbatch{id}.spm.spatial.coreg.estwrite.roptions.interp = 4; % 0: no interpolation 4: 4th degree B-spline interpolation
            matlabbatch{id}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
            matlabbatch{id}.spm.spatial.coreg.estwrite.roptions.mask = 0;
            matlabbatch{id}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
            spm_jobman('run', matlabbatch(id));
        end
        
        %% Image calculation (e.g. merging)
        if ismember('merge', job(j).operation)
            % ROI folder
            if strfind(job(j).ROI{1}, 'w')
                ROIfolder = fullfile(get_path('project'), get_folder(subject, 'r'), ...
                    'fMRI', 'ROI', job(j).atlas); % 'coregistered'
            else
                ROIfolder = fullfile(get_path('project'), get_folder(subject, 'r'), ...
                    'fMRI', 'ROI', job(j).atlas, 'coregistered'); 
            end
            for r=1:nROIs
                file = {}; % size not known beforehand...
                for m=1:2 % files to be merged
                    file(m,:) = getfname(ROIfolder, job(j).ROI{r,m})';
                end
                
                % Merge ROIs
                for i=1:size(file, 2)
                    id = length(matlabbatch) + 1;
                    for hemi=1:2
                        % work-around to avoid NaN whilst merging masks
                        h = spm_vol(fullfile(ROIfolder, file{hemi,i}));
                        V = spm_read_vols(h);
                        V(isnan(V)) = 0;
                        spm_write_image(h, V, fullfile(ROIfolder, file{hemi,i}));
                    end
                    matlabbatch{id}.spm.util.imcalc.input = fullfile(ROIfolder, file(:,i));
                    [pathstr, filename, ext] = fileparts(file{1,i});
                    str = strsplit(job(j).ROI{r,3}, '-');
                    if strfind(job(j).ROI{1}, 'w')
                        matlabbatch{id}.spm.util.imcalc.output = [strrep(filename, ...
                            str{1}, str{2}) '.img']; % output file with pattern replacement nii/img
                    else
                        matlabbatch{id}.spm.util.imcalc.output = [strrep(filename, ...
                            str{1}, str{2}) '.nii']; % output file with pattern replacement nii/img
                    end
                    matlabbatch{id}.spm.util.imcalc.outdir{1} = ROIfolder;
                    matlabbatch{id}.spm.util.imcalc.expression = 'i1 | i2'; % merging
                    matlabbatch{id}.spm.util.imcalc.options.dmtx = 0;
                    matlabbatch{id}.spm.util.imcalc.options.mask = 0;
                    matlabbatch{id}.spm.util.imcalc.options.interp = 1;
                    matlabbatch{id}.spm.util.imcalc.options.dtype = 4;
                    spm_jobman('run', matlabbatch(id));
                end
            end
        end
    end
end