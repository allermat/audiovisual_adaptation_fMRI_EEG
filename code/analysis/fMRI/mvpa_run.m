function mvpa_run

addpath(fullfile(get_path('toolbox'), 'spm12'));

% Load subject info
subjects = subject_info;
subjects = subjects([1:5]);

% Session info for each subject
session = {{1:4 [1 4] 2:3 1 2 3 4 1:2 3:4}; % subj-1
    {6:9 [6 9] 7:8 6 7 8 9 6:7 8:9}; ... % subj-2
    {6:9 [6 9] 7:8 6 7 8 9 6:7 8:9}; ... % subj-3
    {6:9 [6 9] 7:8 6 7 8 9 8:9}; ... % subj-4
    {6:9 [6 9] 7:8 6 7 8 9 6:7 8:9}}; % subj-5

% Set general path
set_path;

for ss=1:length(subjects)
    
    subject = subjects(ss);
    
    for s=1:length(session{ss})
        
        se = find(ismember(subjects(ss).fMRI.session, session{ss}{s}));
        subject.fMRI.session = session{ss}{s};
        
        % -------- Specify decoding configurations and load defaults ---------
        
        % Extra defined fields
        cfg = struct();
        cfg.subject = subject;
        cfg.decoding.machine = 'nu-SVR';
        cfg.learning_curve = 0;
        cfg.example_selection = {};
        % CV and generalization scheme
        % 1 [dirname=resp]: train on pretest (resp+no resp), test on pretest (resp+no resp) & posttest (no resp) 
        % 2 [dirname=mixresp]: train on pretest (resp+no resp), test on pretest (no resp) & posttest (no resp) 
        % 3 [dirname=noresp]: train on pretest (no resp), test on pretest (no resp) & posttest (no resp)
        % 4 [dirname=respresp]: train on pretest (resp+no resp), test on pretest (resp+no resp) & posttest (resp+no resp)
        % 5 [dirname=responly]: train on pretest (resp), test on pretest (resp) & posttest (resp)
        cfg.design.scheme = 3; 
        dirname = {'resp' 'mixresp' 'noresp' 'respresp' 'responly'};
        
        % Toolbox default fields
        cfg.analysis = 'roi';
        cfg.roi.fileindex = 1:5;
        cfg.results.overwrite = 1;
        cfg.results.filestart = [cfg.decoding.machine '_ROI_5_se' sprintf('_%d', se)];
        cfg.results.output = {'predicted_labels'};
        cfg.spmsubdir = 'run_merged2';
        cfg.results.dir = fullfile(get_path('project'), ...
            get_folder(cfg.subject, 'r'), 'fMRI', 'mvpa', 'run_merged2', ...
            [dirname{cfg.design.scheme}, '_scale_multivar_Euc', '_4_fold_cv']);
        mkdir(cfg.results.dir);
        cfg.basic_checks.DoubleFilenameEntriesOk = 1;
        cfg.plot_selected_voxels = 500;
        cfg.plot_design = 1;
        
        % Load defaults
        cfg = decoding_defaults(cfg);
        cfg.parameter_selection.decoding = cfg.decoding;
        
        % ------------------- Run decoding ---------------------
        
        % Read SPM file
        LUT = mvpa_SPM_read(cfg);
        
        % Read files in
        cfg.files = mvpa_make_files(cfg, LUT);
        
        % Calculate design
        cfg.design = make_design_cv_misc(cfg);
        
        % Load mask files (beta images are no longer loaded here)
        mask = mvpa_feat_read(cfg);
        
        % Run decoding (beta images are loaded here whilst also applying
        % multivariate noise normalization)
        mvpa_decode_SVM(mask, cfg);
    end
end