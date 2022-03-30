function w = mvpa_decode_SVM(mask, cfg, fid)

% Load SPM for mutivariate normalization
spmdir = fullfile(get_path('project'), get_folder(cfg.subject, 'r'), ...
    'fMRI', '1st level', cfg.spmsubdir);
load(fullfile(spmdir, 'SPM.mat'));

for rr=1:numel(cfg.roi.mask.fname)
    % Multivariate noise normalization cf. Diedrichsen 2016
    if strfind(cfg.results.dir, 'multivar')
        fname = strrep(mask.ROI(rr).hdr.fname, '.nii', '_m2u_hat.mat');
        if exist(fname, 'file')
            load(fname);
        else
            fprintf('multivariate normalization...');
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
                xY(i,:) = spm_sample_vol(VolIn(i), mask.ROI(rr).XYZ(1,:), ...
                    mask.ROI(rr).XYZ(2,:), mask.ROI(rr).XYZ(3,:), 0);
            end
            % Estimate beta images from raw data whilst also applying
            %  multivariate noise normalization
            [u_hat,resMS,Sw_hat,beta_hat] = rsa.spm.noiseNormalizeBeta(xY, SPM);
            save(fname, 'u_hat');
        end
        data = [];
        for i=1:size(cfg.files.LUT.betaid, 1)
            data(end+1,:) = mean(u_hat(cfg.files.LUT.betaid(i,:),:), 1); % mean from same run estimates
        end
        fprintf('done\n');
        clear xY;
    else
        data = [];
        for i=1:size(cfg.files.LUT.rowid, 1)
            data(end+1,:) = mean(feat{rr}(cfg.files.LUT.rowid(i,:),:), 1); % mean from same run estimates
        end
    end

    % Feature selection preallocation
    if exist('fid', 'var')
        featid = fid{rr};
    else
        featid = true(1, size(data, 2));
    end
    
    % Feature scaling for all data or separately for chunks
    switch cfg.scale.estimation
        case 'all'
            data = mvpa_norm_calc(data, cfg);
        case 'separate' % remove run specific mean condition, cf. Diedrichsen 2016, 2017
            chunks = unique(cfg.files.LUT.overall_run);
            for i=1:numel(chunks)
                data(cfg.files.chunk==(chunks(i)),:) = ...
                    mvpa_norm_calc(data(cfg.files.chunk==(chunks(i)),:), cfg);
            end
        case 'block'
            chunks = unique(cfg.files.LUT.blocktype);
            for i=1:numel(chunks)
                data(ismember(cfg.files.LUT.blocktype, (chunks{i})),:) = ...
                    mvpa_norm_calc(data(ismember(cfg.files.LUT.blocktype, ...
                    (chunks{i})),:), cfg);
            end
    end
    
    % Image scaling
    switch cfg.scale.image
        case 'eucledian' % recommended by PRONTO
            for i=1:size(data, 1)
                data(i,:) = data(i,:) / norm(data(i,:));
            end
        case 'z'
            data = zscore(data, 0, 2);
        case 'mean'
            data = data - repmat(mean(data, 2), 1, size(data, 2));
    end
    
    % Initialize results
    results.predicted_labels.output(rr).predicted_labels = cell(1, size(cfg.design.train, 2));
    
    % Cross-validation loop
    cfg.design.train = logical(cfg.design.train);
    cfg.design.test = logical(cfg.design.test);
    for i=1:size(cfg.design.train, 2)
        
        idtr = find(cfg.design.train(:,i) > 0);
        idte = find(cfg.design.test(:,i) > 0);
        
        % Training and test data
        trdata = data(idtr,:);
        tedata = data(idte,:);
        
        % Traning and test labels
        trlabel = cfg.design.label(idtr);
        telabel = cfg.design.label(idte);
        
        % Feature scaling independently for training and test data
        if strcmp(cfg.scale.estimation, 'across')
            [trdata, scaleparams] = mvpa_norm_calc(trdata, cfg);
            tedata = mvpa_norm_calc(tedata, cfg, scaleparams);
        end
                
        % Train machine
        parameters = cfg.decoding.train.(cfg.decoding.method).model_parameters;
        model = svmtrain(trlabel, trdata(:,featid), parameters);
        
        % Make predictions
        [results.predicted_labels.output(rr).predicted_labels{i}, acc, ...
            decvalue] = svmpredict(telabel, tedata(:,featid),  model);
        w{rr}(:,i) = model.SVs' * model.sv_coef;    
    end
    clear featid;
end

save(fullfile(cfg.results.dir, [cfg.results.filestart ...
    '_predicted_labels.mat']), 'results');
save(fullfile(cfg.results.dir, [cfg.results.filestart '_cfg.mat']), 'cfg');

% LIBSVM options:
% -s svm_type : set type of SVM (default 0)
% 	0 -- C-SVC
% 	1 -- nu-SVC
% 	2 -- one-class SVM
% 	3 -- epsilon-SVR
% 	4 -- nu-SVR
% -t kernel_type : set type of kernel function (default 2)
% 	0 -- linear: u'*v
% 	1 -- polynomial: (gamma*u'*v + coef0)^degree
% 	2 -- radial basis function: exp(-gamma*|u-v|^2)
% 	3 -- sigmoid: tanh(gamma*u'*v + coef0)
% -d degree : set degree in kernel function (default 3)
% -g gamma : set gamma in kernel function (default 1/num_features)
% -r coef0 : set coef0 in kernel function (default 0)
% -c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)
% -n nu : set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)  can only vary on (0,1] cf. Schoelkopf
% -p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)
% -m cachesize : set cache memory size in MB (default 100)
% -e epsilon : set tolerance of termination criterion (default 0.001)
% -h shrinking: whether to use the shrinking heuristics, 0 or 1 (default 1)
% -b probability_estimates: whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 0)
% -wi weight: set the parameter C of class i to weight*C, for C-SVC (default 1)
%
% The k in the -g option means the number of attributes in the input data.
% -v option does cross validation