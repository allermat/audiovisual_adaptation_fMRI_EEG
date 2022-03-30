function out = run_mvpa_fmri

% Define specific analysis
subfolder = fullfile('run_merged2', 'noresp_scale_multivar_Euc_4_fold_cv');
filestart = 'nu-SVR_ROI_5_se_1_2_3_4';

% Load subject info
subjects = subject_info;
subjects = subjects([1:5]);

for ss=1:numel(subjects)
    % Load data and configuration file
    datafolder = fullfile(get_path('project'), get_folder(subjects(ss), 'r'), 'fMRI', 'mvpa', subfolder);
    load(fullfile(datafolder, [filestart '_predicted_labels.mat']));
    load(fullfile(datafolder, [filestart '_cfg.mat']));
    
    nROIs = numel(cfg.roi.mask.name);
    
    cfg.design.test = logical(cfg.design.test);
    
    for i=1:nROIs
        maxset = max(cfg.design.set);
        maxcol = find(cfg.design.set == maxset);
        maxlabel{i} = arrayfun(@(x) cfg.design.label(cfg.design.test(:,x)>0,x), maxcol, 'un', 0);
        maxpredlabel{i} = results.predicted_labels.output(i).predicted_labels(maxcol);
        blocktype{i} = arrayfun(@(x) cfg.files.blocktype(cfg.design.test(:,x)>0), maxcol, 'un', 0);
    end
    
    for i=1:nROIs
        maxlabel{i} = arrayfun(@(x) cell2mat(maxlabel{i}(:,x)), 1:size(maxlabel{i}, 2), 'un', 0);
        maxpredlabel{i} = arrayfun(@(x) cell2mat(maxpredlabel{i}(:,x)), 1:size(maxpredlabel{i}, 2), 'un', 0);
        blocktype{i} = arrayfun(@(x) cat(1, blocktype{i}{:,x}), 1:size(blocktype{i}, 2), 'un', 0);
        
        cv_label = cat(2, maxlabel{:,i});
        cv_predlabel = cat(2, maxpredlabel{:,i});
        cv_blocktype = cat(2, blocktype{:,i});
        
        % Take average CV result if same examples were generalized multiple times!
        post_label = maxlabel{i}{1};
        post_predlabel = mean(cat(2, cv_predlabel{:}), 2);
        post_blocktype = blocktype{i}{1};
        
        % Calculate spatial encoding index of pretest data
        out.data(ss).accuracy(i,1) = mean(atanh(cellfun(@(x,y,z) ...
            corr(x(ismember(z, 'pretest')),y(ismember(z, 'pretest'))), ...
            cv_label, cv_predlabel, cv_blocktype))); % tanh for corr!
               
        % Calculate recalibration index of posttest data
        StimLevels = unique(cat(1, cv_label{:}));
        adapt_all = {'posttest-radapt' 'posttest-ladapt'};
        for sl=1:numel(StimLevels)
            for ad=1:numel(adapt_all)
                id = ismember(post_blocktype, adapt_all{ad}) ...
                    & post_label == StimLevels(sl);
                out.data(ss).OutOfNum{i,ad+1}(sl,1) = sum(id);
                out.data(ss).NumPos{i,ad+1}(sl,1) = sum(post_predlabel(id) > 0);
            end
            out.data(ss).recal_index(i,sl) = out.data(ss).NumPos{i,2}(sl) ...
                / out.data(ss).OutOfNum{i,2}(sl) - ...
                out.data(ss).NumPos{i,3}(sl) / out.data(ss).OutOfNum{i,3}(sl);
        end
        
        % Concatenate data
        all_label = cat(1, cv_label{:});
        all_predlabel = cat(1, cv_predlabel{:});
        all_blocktype = cat(1, cv_blocktype{:});
        preid = ismember(all_blocktype, 'pretest');
        label = [all_label(preid) all_predlabel(preid)];
        
        % Calculate number of positive and out of number responses
        for sl=1:numel(StimLevels)
            out.data(ss).OutOfNum{i,1}(sl,1) = sum(all_label(preid)==StimLevels(sl));
            out.data(ss).NumPos{i,1}(sl,1) = sum(all_predlabel(preid ...
                & all_label==StimLevels(sl)) > 0);
        end
    end
end

% Save ROIs into output
out.ROI = cfg.roi.mask.name;
out.ROI{ismember(out.ROI, 'A1')} = 'HG';
out.ROI{ismember(out.ROI, 'IPS0-5')} = 'IPS';