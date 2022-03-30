function [gimg,aloc,conds,roiNames] = lme_getData(varargin)
% Load and prepare fMRI data for BOLD and LME analysis. 
% 
% Input: 
%   hemisphere: from which hemisphere the data should be loaded ('lh','rh')
%       If 'mean_lhrh', then the data are averaged across hemispheres such that
%       the order of conditions are flipped in the right hemisphere. 
% Output: 
%   gimg: 4D array of BOLD acrivity: 
%          nlocs (-12 -5 -2 0 2 5 12) x 
%          nconds (pre, postAV, postVA) x 
%          nsubjects x 
%          nROIs
%   aloc: auditory stimulus locations (-12 -5 -2 0 2 5 12)
%   conds: conditions (pre, postAV, postVA)
%   roiNames: names of ROIs (HG, hA, IPL, IPS, FEF)
%

% Parsing input
p = inputParser;
validHemispheres = {'mean_lhrh','lh','rh'};
addParameter(p,'hemisphere','mean_lhrh',@(x) ismember(x,validHemispheres));
parse(p,varargin{:});
hemisphere = p.Results.hemisphere;

% Load subject info
subjects = subject_info;
subjects = subjects([1:5]);
nsubjects = length(subjects);
nROIs = 10;

filestart = 'con';
fileend = '';
aloc = [-12 -5 -2 0 2 5 12];
nlocs = numel(aloc);
nconds = 3;

gimg = zeros(nconds, nlocs, nsubjects, nROIs);

for ss=1:nsubjects
    fprintf('\n%s\n', subjects(ss).id);
    
    % Configuration
    cfg.subject = subjects(ss);
    cfg.spmsubdir = 'run_merged2';
    cfg.roi.fileindex = 6:15; % left (6:10) and right (11:15) hemispheres
    cfg = decoding_defaults(cfg);
    
    % SPM output folder
    spmdir = fullfile(get_path('project'), get_folder(subjects(ss), 'r'), 'fMRI', ...
        '1st level', cfg.spmsubdir);
    
    % Get 'no-response' files and corresponding T contrast from SPM
    cfg.files.name = arrayfun(@(x) fullfile(spmdir, getfname(spmdir, ...
        sprintf('%s_%04d%s.nii', filestart, x, fileend))), 2:2:42);
    cfg.files.name{end+1} = fullfile(spmdir, 'spmT_0044.nii');
    
    % Get mask file
    cfg.files.mask =  cellfun(@(x,y) fullfile(get_path('project'), ...
        get_folder(subjects(ss), 'r'), 'fMRI', 'ROI', x, 'coregistered', ...
        [y cfg.roi.extension]), cfg.roi.mask.atlas, cfg.roi.mask.fname, 'un', 0);
    
    % Read in masked images
    [mask, img] = mvpa_feat_read(cfg);
    
    % Select top20 voxels
    [~, sortid] = cellfun(@(x) sort(x(end,:), 2, 'descend'), img, 'un', 0);
    sortid = cellfun(@(x) x(1:20), sortid, 'un', 0);
    for i=1:nROIs
        img{i} = img{i}(:,sortid{i});
        img{i}(end,:) = [];
    end
    
    % Reshape data to pre/postAV/postVA
    gimg(:,:,ss,:) = reshape(cell2mat(cellfun(@(x) reshape(nanmean(x, 2), ...
        nconds, nlocs, []), img, 'un', 0)), nconds, nlocs, nROIs); 
end

switch hemisphere
    case 'mean_lhrh'
         % flip direction
        [gimg(2,:,:,6:10), gimg(3,:,:,6:10), gimg(1,:,:,6:10)] = deal(...
            gimg(3,7:-1:1,:,6:10), gimg(2,7:-1:1,:,6:10), gimg(1,7:-1:1,:,6:10));
        gimg = reshape(gimg, nconds, nlocs, nsubjects, 5, 2);
        gimg = mean(gimg, 5);
    case 'lh'
        gimg = gimg(:,:,:,1:5);
    case 'rh'
        gimg = gimg(:,:,:,6:10);
end

roiNames = {'HG','hA','IPL','IPS','FEF'};
conds = {'pre','postAV','postVA'};
gimg = permute(gimg, [2 1 3 4]);


end