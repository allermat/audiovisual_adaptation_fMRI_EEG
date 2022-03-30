function [pcmDataAll,LUTAll] = pcm_getData(varargin) 

validDistFun = {'Euclidean','Mahalanobis','Correlation','cvEuclidean',...
    'cvMahalanobis'};
validNormMode = {'runwise','overall'};
validBlockTypes = {'pretest','posttest-ladapt' 'posttest-radapt'};
validOutput = {'mean','byChunk','bySession'};


p = inputParser;
addParameter(p,'distfun','Mahalanobis',@(x) ismember(x,validDistFun));
addParameter(p,'blocktype',{{'pretest'}},@(x) all(cellfun(@(y) all(ismember(y,validBlockTypes)),x)));
addParameter(p,'ismultinorm',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'normmode','overall',@(x) ismember(x,validNormMode));
addParameter(p,'iseucnorm',false,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'byHemisphere',false,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'output','mean',@(x) ismember(x,validOutput));
addParameter(p,'topNvoxelTval',[],@(x) validateattributes(x,{'numeric'},...
    {'scalar','positive'}));
addParameter(p,'topNvoxelSVM',[],@(x) validateattributes(x,{'numeric'},...
    {'scalar','positive'}));

parse(p,varargin{:})

distfun = p.Results.distfun;
ismultinorm = p.Results.ismultinorm;
normmode = p.Results.normmode;
iseucnorm = p.Results.iseucnorm;
output = p.Results.output;
blocktype = p.Results.blocktype;
byHemisphere = p.Results.byHemisphere;
topNvoxelTval = p.Results.topNvoxelTval;
topNvoxelSVM = p.Results.topNvoxelSVM;

% if strfind(distfun, 'Mahalanobis')
%     ismultinorm = 1;
% end

% Load subject info
subjects = subject_info;
subjects = subjects([1:5]);
nsubjects = numel(subjects);

if byHemisphere
    iROI = [1:10]; % 1:5
else
    iROI = [1:5]; % 1:5
end

if ~isempty(topNvoxelSVM)
    featid = load(fullfile(get_path('project'),'data','group','fMRI','RSA',...
        'featid.mat'),sprintf('featid_%d',topNvoxelSVM));
    featid = featid.(sprintf('featid_%d',topNvoxelSVM));
end
nROIs = numel(iROI);
allpath = {};
for bb=numel(blocktype)
    
    beta = cell(nsubjects, nROIs);
    [pcmDataAll,LUTAll] = deal(cell(nsubjects, nROIs));
    subjectRDMs = cell(nROIs, 1);
    
    for ss=1:nsubjects
        subfolder = fullfile('run_merged2', 'noresp_scale_multivar_Euc_4_fold_cv');
        if byHemisphere
            filestart = ['nu-SVR_ROI_10_se' sprintf('_%d', 1:4)];
        else
            filestart = ['nu-SVR_ROI_5_se' sprintf('_%d', 1:4)];
        end
        datafolder = fullfile(get_path('project'), get_folder(subjects(ss), 'r'), ...
            'fMRI', 'mvpa', subfolder);
        load(fullfile(datafolder, [filestart '_cfg.mat']));
        for i=1:nROIs
            % Update filename for cross-platform compatibility
            cfg.files.mask{i} = strrep(cfg.files.mask{i},'D:\','d:\'); 
            cfg.files.mask{i} = strrep(cfg.files.mask{i}, ...
                'd:\A_Recalibration_fMRI', get_path('project'));
        end
        allpath = [allpath; fullfile(datafolder, [filestart '_cfg.mat']); cfg.files.mask];

        if ~isempty(topNvoxelTval)
            % Configuration
            cfg2.subject = subjects(ss);
            cfg2.spmsubdir = 'run_merged2';
            if byHemisphere
                cfg2.roi.fileindex = 6:15;
            else
                cfg2.roi.fileindex = 1:5;
            end
            cfg2 = decoding_defaults(cfg2);
            spmdir = fullfile(get_path('project'),get_folder(subjects(ss),'r'),...
                'fMRI','1st level', cfg2.spmsubdir);
            cfg2.files.name{1} = fullfile(spmdir, 'spmT_0044.nii'); % noresp
            cfg2.files.mask =  cellfun(@(x,y) fullfile(get_path('project'), ...
                get_folder(subjects(ss), 'r'), 'fMRI', 'ROI', x, 'coregistered', ...
                [y cfg2.roi.extension]), cfg2.roi.mask.atlas, cfg2.roi.mask.fname, 'un', 0);
            % Read in masked images
            [mask, img] = mvpa_feat_read(cfg2);
            % Sort beta image and choose the top x voxel
            [~, sortid] = cellfun(@(x) sort(x(end,:), 2, 'descend'), img, 'un', 0);
            sortid = cellfun(@(x) x(1:topNvoxelTval), sortid, 'un', 0);
        elseif ~isempty(topNvoxelSVM)
            % Sort SVM weight image and choose the top x voxel
            [~, sortid] = cellfun(@(x) sort(x(end,:), 2, 'descend'), featid(ss,:), 'un', 0);
            sortid = cellfun(@(x) x(1:topNvoxelSVM), sortid, 'un', 0);
        end
        
        % Read in beta
        for rr=1:nROIs
            if strcmp(cfg.spmsubdir, 'run_merged2') && ismultinorm
                fname = strrep(cfg.files.mask{iROI(rr)}, '.nii', '_m2u_hat');
            else
                error('This setting is no longer supported.');
            end
            if strcmp(normmode, 'runwise')
               fname = [fname '_run']; 
            end
            if exist([fname '.mat'], 'file')
                allpath = [allpath; [fname '.mat']];
                load([fname '.mat'])
            else
                % Multivariate noise normalization cf. Diedrichsen 2016
                fprintf('multivariate normalization...');
                spmdir = fullfile(get_path('project'), get_folder(cfg.subject, 'r'), ...
                    'fMRI', '1st level', cfg.spmsubdir);
                load(fullfile(spmdir, 'SPM.mat'));
                VolIn = SPM.xY.VY;
                for i = 1:numel(VolIn)
                    % Backward compatibility: propagate scaling (see spm_fmri_spm_ui.m)
                    VolIn(i).private.dat.scl_slope = VolIn(i).pinfo(1);
                    VolIn(i).private.dat.scl_inter = VolIn(i).pinfo(2);
                end
                for i=1:length(VolIn)
                    % Update filename for cross-platform compatibility
                    VolIn(i).fname = strrep(VolIn(i).fname, ...
                        'd:\A_Recalibration_fMRI', get_path('project'));
                    xY(i,:) = spm_sample_vol(VolIn(i), mask.ROI(iROI(rr)).XYZ(1,:), ...
                        mask.ROI(iROI(rr)).XYZ(2,:), mask.ROI(iROI(rr)).XYZ(3,:), 0);
                end
                % Estimate beta images from raw data whilst also applying
                % multivariate noise normalization
                [u_hat,resMS,Sw_hat,beta_hat] = rsa.spm.noiseNormalizeBeta(xY, ...
                    SPM, [], 'normmode', normmode);
                save(fname, 'u_hat');
                clear xY;
                fprintf('done\n');
            end
            
            for i=1:size(cfg.files.LUT.betaid, 1)
                % Mean from same run estimates
                beta{ss,rr}(end+1,:) = mean(u_hat(cfg.files.LUT.betaid(i,:),:), 1);
            end
            
            % Eucledian normalization
            if iseucnorm
                for i=1:size(beta{ss,rr}, 1)
                    beta{ss,rr}(i,:) = beta{ss,rr}(i,:) / norm(beta{ss,rr}(i,:));
                end
            end
            
            % Calculate RDM
            aloc = unique(cfg.files.LUT.aloc);
            nblocks = numel(blocktype{bb});
            if strfind(distfun, 'cv')
                clear run;
                for b=1:nblocks
                    run{b} = unique(cfg.files.LUT.overall_run(ismember(cfg.files.LUT.blocktype, ...
                        blocktype{bb}{b})));
                end
                nexamples = size(beta{ss,rr}, 1);
                [partition, conditionVec] = deal(zeros(nexamples, 1));
                for i=1:nexamples
                    if cfg.files.LUT.resp(i) == 1 && ismember(cfg.files.LUT.blocktype{i}, blocktype{bb})
                        switch cfg.files.LUT.blocktype{i}
                            case {'posttest-ladapt' 'pretest'}
                                conditionVec(i) = find(aloc == cfg.files.LUT.aloc(i));
                            case 'posttest-radapt'
                                conditionVec(i) = find(aloc == cfg.files.LUT.aloc(i)) + 7;
                        end
                        partition(i) = find(run{ismember(blocktype{bb}, ...
                            cfg.files.LUT.blocktype{i})} == cfg.files.LUT.overall_run(i));
                    end
                end
                subjectRDMs{rr}(:,:,ss) = squareform(rsa.distanceLDC(beta{ss,rr}, partition, conditionVec));
            else
                if strcmp(output,'byChunk')
                    pcm_data = cell(nblocks, 1);
                    nChunks = length(unique(cfg.files.LUT.chunk(~isnan(cfg.files.LUT.chunk))));
                    for b=1:nblocks
                        for ch = 1:nChunks
                            for al=1:numel(aloc)
                                pcm_data{b}(al,ch,:) = mean(...
                                    beta{ss,rr}(ismember(cfg.files.LUT.blocktype,blocktype{bb}{b}) & ...
                                    cfg.files.LUT.resp == 0 & ...
                                    cfg.files.LUT.aloc == aloc(al) & ...
                                    cfg.files.LUT.chunk == ch,:), 1);
                            end
                        end
                    end
                    pcm_data = cat(1, pcm_data{:});
                    pcm_data = reshape(pcm_data,size(pcm_data,1)*size(pcm_data,2),[]);
                    pcmDataAll{ss,rr} = pcm_data;
                elseif strcmp(output,'bySession')
                    pcm_data = cell(nblocks,1);
                    LUT = repmat({table()},nblocks,1);
                    sessions = unique(cfg.files.LUT.spm_session);
                    for b=1:nblocks
                        for sn = 1:numel(sessions)
                            for al=1:numel(aloc)
                                isSel = ismember(cfg.files.LUT.blocktype,blocktype{bb}{b}) & ...
                                    cfg.files.LUT.resp == 0 & ...
                                    cfg.files.LUT.aloc == aloc(al) & ...
                                    cfg.files.LUT.spm_session == sessions(sn);
                                pcm_data{b} = cat(1,pcm_data{b},beta{ss,rr}(isSel,:));
                                LUT{b} = cat(1,LUT{b},cfg.files.LUT(isSel,:));
                            end
                        end
                    end
                    temp = cat(1,pcm_data{:});
                    if ~isempty(topNvoxelTval) || ~isempty(topNvoxelSVM)
                        temp = temp(:,sortid{rr});
                    end
                    pcmDataAll{ss,rr} = temp;
                    LUTAll{ss,rr} = cat(1,LUT{:});
                elseif strcmp(output,'mean')
                    pcm_data = cell(nblocks, 1);
                    for b=1:nblocks
                        for al=1:numel(aloc)
                            pcm_data{b}(al,:) = mean(...
                                beta{ss,rr}(ismember(cfg.files.LUT.blocktype,blocktype{bb}{b}) & ...
                                cfg.files.LUT.resp == 0 & ...
                                cfg.files.LUT.aloc == aloc(al),:), 1);
                        end
                    end
                    pcm_data = cat(1, pcm_data{:});
                    pcmDataAll{ss,rr} = pcm_data;
                end
            end
        end
    end
end