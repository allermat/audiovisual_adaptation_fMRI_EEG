function [mask, feat] = mvpa_feat_read(cfg)

% --------- Operations with whole brain global mask ----------

spmdir = fullfile(get_path('project'), get_folder(cfg.subject, 'r'), ...
    'fMRI', '1st level', cfg.spmsubdir);
mask.global.hdr = spm_vol(char(fullfile(spmdir, 'mask.nii')));
mask.global.img = spm_read_vols(mask.global.hdr);

% ----------- Operations with ROI masks -----------

% ROI headers
for i=1:length(cfg.files.mask)
    mask.ROI(i).hdr = spm_vol(char(cfg.files.mask{i}));
end

% Check whether global and ROI masks are in same space
hdr = cat(1, mask.ROI.hdr);
sts = spm_check_orientations([mask.global.hdr; hdr]);
if sts ~= 1
    error('Images not in same space!');
end

% Create mask in XZY format (both world and voxel coordinates)
linid = find(mask.global.img > 0);
[X, Y, Z] = ind2sub(size(mask.global.img), linid);
mask.global.XYZ = [X'; Y'; Z']; % XYZ format
mask.global.size = size(mask.global.XYZ, 2);
mask.global.XYZmm = mask.global.hdr.mat(1:3,:) * [mask.global.XYZ; ...
    ones(1, mask.global.size)]; % voxel to world transformation
% Combine masks
for i=1:length(mask.ROI)
    xY.def = 'mask';
    xY.spec = char(cfg.files.mask{i});
    [xY, mask.ROI(i).XYZmm, j] = spm_ROI(xY, mask.global.XYZmm);
    mask.ROI(i).XYZ = mask.global.XYZ(:,j);
    mask.ROI(i).size = size(mask.ROI(i).XYZ, 2);
end

if nargout == 2
    % ------------- Create feature set ---------------
    
    % Mask each image by each ROI and create overall feature set (images x voxel)
    fprintf('preparing feature set...');
    [C, ia, ic] = unique(cfg.files.name); % get only non-repetitive file names
    V = spm_vol(char(C));
    feat = cell(1, length(mask.ROI));
    for i=1:length(mask.ROI)
        data = spm_get_data(V, mask.ROI(i).XYZ);
        feat{i} = data(ic,:); % add repetitions back
        clear spm_sample_vol % free up memory
    end
    fprintf('done\n');
end