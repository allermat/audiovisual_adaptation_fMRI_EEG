function RDM = get_RDMs_fMRI(varargin)

p = inputParser();

validDistFuns = {'Euclidean','Mahalanobis','Correlation','cvEuclidean','cvMahalanobis'};
validRdmLayouts = {'pre','post','prepost'};

addParameter(p,'distfun','Mahalanobis',@(x) ismember(x,validDistFuns));
addParameter(p,'rdmLayout',{'pre','post'},@(x) all(ismember(x,validRdmLayouts)));
addParameter(p,'normmode','overall',@(x) ismember(x,{'runwise','overall'}));
addParameter(p,'ismultinorm',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'iseucnorm',true,@(x) validateattributes(x,{'logical'},{'scalar'}));

parse(p,varargin{:});

distfun = p.Results.distfun;
rdmLayout = p.Results.rdmLayout;
normmode = p.Results.normmode;
ismultinorm = p.Results.ismultinorm;
iseucnorm = p.Results.iseucnorm;

saveSubjectRDMs = false;

addpath(genpath(fullfile(get_path('toolbox'), 'rsatoolbox')));

% Load subject info
subjects = subject_info;
subjects = subjects([1:5]);
nsubjects = numel(subjects);

iROI = [1:5]; % 1:5
nROIs = numel(iROI);

allpath = {};
blocktype = cell(size(rdmLayout));
for i = 1:numel(blocktype)
    switch rdmLayout{i}
        case 'pre'
            blocktype{i} = {'pretest'};
        case 'post'
            blocktype{i} = {'posttest-ladapt','posttest-radapt'};
        case 'prepost'
            blocktype{i} = {'pretest','posttest-ladapt' 'posttest-radapt'};
    end
end

for bb=1:numel(blocktype)
    
    beta = cell(nsubjects, nROIs);
    subjectRDMs = cell(nROIs, 1);
    
    for ss=1:nsubjects
        subfolder = fullfile('run_merged2', 'noresp_scale_multivar_Euc_4_fold_cv');
        filestart = ['nu-SVR_ROI_5_se' sprintf('_%d', 1:4)]; % 6 IC2n
        
        datafolder = fullfile(get_path('project'), get_folder(subjects(ss), 'r'), 'fMRI', 'mvpa', subfolder);
        load(fullfile(datafolder, [filestart '_cfg.mat']));
        for i=1:nROIs
            cfg.files.mask{i} = strrep(cfg.files.mask{i}, ... % update path...
                'd:\A_Recalibration_fMRI', get_path('project'));
        end
        allpath = [allpath; fullfile(datafolder, [filestart '_cfg.mat']); cfg.files.mask];

        % Load ROI mask
%         mask = mvpa_ROI_read(cfg);
        
        % Read in beta
        for rr=1:nROIs
            if strcmp(distfun, 'Mahalanobis') && ismultinorm == 0
                fname = strrep(cfg.files.mask{iROI(rr)}, '.nii', '_Sw_reg');
            elseif strcmp(cfg.spmsubdir, 'run_merged2') && ismultinorm
                fname = strrep(cfg.files.mask{iROI(rr)}, '.nii', '_m2u_hat');
            elseif ismultinorm
                fname = strrep(cfg.files.mask{iROI(rr)}, '.nii', '_mu_hat'); % _u_hat_run _beta_hat
            else
                fname = strrep(cfg.files.mask{iROI(rr)}, '.nii', '_mbeta_hat');
            end
            if strcmp(normmode, 'runwise')
               fname = [fname '_run']; 
            end
            if exist([fname '.mat'], 'file')
                allpath = [allpath; [fname '.mat']];
                load(fname)
                if strcmp(distfun, 'Mahalanobis') && ismultinorm == 0
                    load(strrep(cfg.files.mask{iROI(rr)}, '.nii', '_mbeta_hat.mat'));
                end
            else
                % Multivariate noise normalization cf. Diedrichsen 2016
                fprintf('multivariate normalization...');
                spmdir = fullfile(get_path('project'), get_folder(cfg.subject, 'r'), 'fMRI', '1st level', cfg.spmsubdir);
                load(fullfile(spmdir, 'SPM.mat'));
                VolIn = SPM.xY.VY;
                for i = 1:numel(VolIn)
                    % Backward compatibility: propagate scaling (see spm_fmri_spm_ui.m)
                    VolIn(i).private.dat.scl_slope = VolIn(i).pinfo(1);
                    VolIn(i).private.dat.scl_inter = VolIn(i).pinfo(2);
                end
                for i=1:length(VolIn)
                    xY(i,:) = spm_sample_vol(VolIn(i), mask.ROI(iROI(rr)).XYZ(1,:), mask.ROI(iROI(rr)).XYZ(2,:), mask.ROI(iROI(rr)).XYZ(3,:), 0);
                end
                if strcmp(distfun, 'Mahalanobis') && ismultinorm == 0
                    [Sw_reg, res, beta_hat] = rsa.get_res_covmatrix(xY, SPM);
                else
                    [u_hat,resMS,Sw_hat,beta_hat] = rsa.spm.noiseNormalizeBeta(xY, SPM, 'normmode', normmode);
                end
                if strcmp(distfun, 'Mahalanobis') && ismultinorm == 0
                    save(fname, 'Sw_reg');
                elseif ismultinorm
                    save(fname, 'u_hat');
                else
                    save(fname, 'beta_hat');
                end
                clear xY;
                fprintf('done\n');
            end
            for i=1:size(cfg.files.LUT.betaid, 1)
                if ismultinorm
                    beta{ss,rr}(end+1,:) = mean(u_hat(cfg.files.LUT.betaid(i,:),:), 1); % mean from same run estimates
                else
                    beta{ss,rr}(end+1,:) = mean(beta_hat(cfg.files.LUT.betaid(i,:),:), 1);
                end
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
                    run{b} = unique(cfg.files.LUT.overall_run(ismember(cfg.files.LUT.blocktype, blocktype{bb}{b})));
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
                        partition(i) = find(run{ismember(blocktype{bb}, cfg.files.LUT.blocktype{i})} == cfg.files.LUT.overall_run(i));
                        %                             partition(i) = cfg.files.LUT.chunk(i);
                    end
                end
                subjectRDMs{rr}(:,:,ss) = squareform(rsa.distanceLDC(beta{ss,rr}, partition, conditionVec));
            else
                rsa_data = cell(nblocks, 1);
                for b=1:nblocks
                    for al=1:numel(aloc)
                        rsa_data{b}(al,:) = mean(beta{ss,rr}(ismember(cfg.files.LUT.blocktype, blocktype{bb}{b}) & cfg.files.LUT.resp==0 & cfg.files.LUT.aloc == aloc(al),:), 1);
                    end
                end
                rsa_data = cat(1, rsa_data{:});
                if strcmp(distfun, 'Mahalanobis') && ismultinorm == 0
                    subjectRDMs{rr}(:,:,ss) = squareform(pdist(rsa_data, lower(distfun), Sw_reg));
                elseif strcmp(distfun, 'Mahalanobis')
                    subjectRDMs{rr}(:,:,ss) = squareform(pdist(rsa_data, 'euclidean')) .^ 2;
                else
                    subjectRDMs{rr}(:,:,ss) = squareform(pdist(rsa_data, lower(distfun)));
                end
            end
        end
    end
    % Save subject RDMs
    if saveSubjectRDMs
        save(sprintf('subjectRDMs_%s_%s',blocktype{bb}{1}(1:3),distfun(1:3)),...
            'subjectRDMs');
    end
    RDM.(rdmLayout{bb}) = subjectRDMs;    
end

end


function RDM = sub_ranktransform(RDM, nblocks)

id = [0 cumsum(repmat(7, 1, nblocks))] + 1;
for k=1:numel(id)-1
    for j=1:numel(id)-1
        RDM(id(k):id(k+1)-1,id(j):id(j+1)-1) = rsa.util.scale01(rsa.util.rankTransform_equalsStayEqual(RDM(id(k):id(k+1)-1,id(j):id(j+1)-1), 1));
    end
end

end