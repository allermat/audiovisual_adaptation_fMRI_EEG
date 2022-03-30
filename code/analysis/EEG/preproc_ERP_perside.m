function preproc_ERP_perside(subID)
%% Parsing input, checking matlab
p = inputParser;

addRequired(p,'subID',@ischar);

parse(p,subID);

subID = p.Results.subID;

%% Preparing file and directory names for the processing pipeline.
expStage = 'final';

% Loading files specifying parameters
trigDef = load(fullfile(get_path('project'),'experiment','trigger_LUT.mat'));
trigDef = trigDef.trigger_LUT;
setupSpec = load(fullfile(get_path('project'),'analysis','EEG','setup_spec.mat'));
setupSpec = setupSpec.setup_spec;
subjectSpec = subject_info;
subjectSpec = subjectSpec(ismember({subjectSpec.id}',subID));
if isempty(subjectSpec)
    error('preproc_MVPA:invalidInput', ...
          'Invalid subject ID');
end
if isempty(fieldnames(subjectSpec.EEG))
    error('preproc_MVPA:missingData', ...
          'Subject specification does not contain EEG sesssions');
end

sourceDir = fullfile(get_path('project'),'data',subID,'EEG','preproc_data','COMM');
destDir = fullfile(get_path('project'),'data',subID,'EEG','preproc_data','ERP');
if ~exist(destDir,'dir');
    mkdir(destDir);
end

dataSession = cell(numel(subjectSpec.EEG.session),1);

fprintf('\n\nProcessing files...\n');

for iSession = 1:numel(subjectSpec.EEG.session)
    
    actSessionFolderName = subjectSpec.session_folder_name{subjectSpec.EEG.session(iSession)};
    rawDataDir = fullfile(get_path('project'),'data',subID,'EEG','raw_data',actSessionFolderName);
    if ~exist(rawDataDir,'dir')
        error('The specified data folder does not exist!');
    end
    behavSourceDir = fullfile(get_path('project'),'data',subID, ...
                              'EEG','behavioural data',actSessionFolderName);
    
    rawDataFileNamesActDay = {subjectSpec.EEG.filespec(iSession).file.name}';
    sourceFileNamesActSession = strcat('fteeg_COMM_',rawDataFileNamesActDay,'.mat');
    artefactFileNamesActDay = strcat('artf_COMM_',rawDataFileNamesActDay,'.mat');
    
    % Cell array for collecting result files
    resultFilesActSession = cell(numel(sourceFileNamesActSession),1);
    
    for iFileActSession = 1:size(sourceFileNamesActSession,1)
        %% Loading the source data
        ftDataReref = load(fullfile(sourceDir,sourceFileNamesActSession{iFileActSession}));
        ftDataReref = ftDataReref.ftDataReref;
        actFileSpec = subjectSpec.EEG.filespec(iSession).file(iFileActSession);
        actFileSpec.adaptdir = subjectSpec.EEG.pretest.match{iSession};
        
        %% low-pass filtering
        cfg = struct();
        cfg.lpfilter = 'yes';
        cfg.lpfreq = 45;
        cfg.lpfiltord = 5;
        cfg.lpfiltdir = 'twopass';
        ftDataLp = ft_preprocessing(cfg,ftDataReref);
        % Clearing unnecessary previous dataset
        ftDataReref = []; 
        
        %% Epoching
        % Trial definition
        cfg = struct();
        cfg.headerfile = fullfile(rawDataDir,[rawDataFileNamesActDay{iFileActSession},'.vhdr']);
        cfg.trialfun = 'ft_trialfun_eventlocked';
        cfg.trialdef = struct();
        cfg.trialdef.blocktype = {'pre-test','post-test'};
        cfg.trialdef.prestim = 0.1;
        cfg.trialdef.poststim = 0.5;
        cfg.trialdef.trigdef = trigDef;
        cfg.trialdef.fileSpec = actFileSpec;
        cfg.trialdef.trig_visonset_corr = setupSpec.trig_visonset_corr_eeg;
        cfg.trialdef.eventtype = 'Stimulus';
        cfg = ft_definetrial(cfg);
        
        % epoching with background correction
        ftDataEp = ft_redefinetrial(cfg,ftDataLp);
        
        % Clearing unnecessary previous dataset
        ftDataLp = []; 
        
        %% Baseline correction
        cfg = struct();
        cfg.demean = 'yes';
        cfg.baselinewindow = [-0.1,0];
        ftDataEp = ft_preprocessing(cfg,ftDataEp);
        
        %% Rejecting trials based on artefacts, eye events and behaviour. 
                
        % Loading artefact data
        dataArtf = load(fullfile(sourceDir,artefactFileNamesActDay{iFileActSession}));
        dataBehav = loadBehavData(behavSourceDir,actFileSpec);
        
        % Window of interest around the onset of the stimulus which has to
        % be free of artefacts and unwanted eye events. 
        winOfInt = [-0.1,0.5];
        
        % Marking bad trials
        % Last column of trialinfo gives the following info:
        % 0 - good
        % 1 - no response
        % 2 - false alarm
        % 3 - wrong hand (not used)
        % 4 - EEG artefact
        % 5 - missing eyetracker data
        % 6 - eyeblink
        % 7 - saccade
        % 8 - wrong fixation location
        badTrials = rejecttrials(ftDataEp,dataBehav,dataArtf,trigDef,winOfInt,behavSourceDir);
        
        ftDataEp.trialinfo = buildtrialinfo([ftDataEp.trialinfo,badTrials],dataBehav,actFileSpec);
        
        resultFilesActSession{iFileActSession} = ftDataEp;
        
        ftDataEp = []; 
        
    end
        
    %% Merging files for the same day if necessary
    if numel(resultFilesActSession) > 1
        ftDataClean = ft_appenddata([],resultFilesActSession{:});
    else
        ftDataClean = resultFilesActSession{1}; %#ok<*NASGU>
    end
    
    resultFilesActSession = [];
    
    dataSession{iSession} = ftDataClean;
    
    ftDataClean = [];
    
end

%% Merging files across sessions if necessary
if numel(dataSession) > 1
    ftDataMerged = ft_appenddata([],dataSession{:});
else
    ftDataMerged = dataSession{1}; %#ok<*NASGU>
end

%% Downsampling
cfg = struct();
cfg.resamplefs = 200;
cfg.detrende = 'no';
ftDataDsampl = ft_resampledata(cfg,ftDataMerged);
% Clearing unnecessary previous dataset
ftDataMerged = [];

%% Averaging over conditions
conditions = trigDef(ismember(trigDef.blocktype,{'pre-test','post-test'}),...
                     {'aloc','blocktype','condition'});
temp = conditions(ismember(conditions.blocktype,'post-test'),:);
conditions.adaptdir = repmat({'L'},size(conditions,1),1);
temp.adaptdir = repmat({'R'},size(temp,1),1);
conditions = cat(1,conditions,temp);
conditions.adaptdir(ismember(conditions.blocktype,'pre-test')) = {''};

condsPerSide = conditions(ismember(conditions.aloc,[-12,-5]),:);
condsPerSide.condition = cell(size(condsPerSide,1),1);
condsPerSide.aloc(condsPerSide.aloc == -12) = -1;
condsPerSide.aloc(condsPerSide.aloc == -5) = 1;
adaptDirs = {'L','R'};
for i = [-1,1]
    condsPerSide.condition(ismember(condsPerSide.blocktype,'pre-test') & ...
        condsPerSide.aloc == i) = ...
        {conditions.condition(sign(conditions.aloc) == sign(i) & ...
        ismember(conditions.blocktype,'pre-test'))};
    for j = 1:2
        condsPerSide.condition(ismember(condsPerSide.blocktype,'post-test') & ...
            condsPerSide.aloc == i) = ...
            {conditions.condition(sign(conditions.aloc) == sign(i) & ...
            ismember(conditions.blocktype,'post-test') & ...
            ismember(conditions.adaptdir,adaptDirs{j}))};
    end
end
% Initializing progress monitor
fprintf('\n\nAveraging tials...\n');
parfor_progress(size(condsPerSide,1));

for iCond = 1:size(condsPerSide,1)
    
    actCond = condsPerSide(iCond,:);
    
    temp = [ftDataDsampl.trialinfo{:}];
    cond = [temp.condition];
    badTrials = [temp.badtrials];
    catchTrials = [temp.catch_trial];
    adaptdir = {temp.adaptdir};
    
    cfg = struct();
    % Selecting good trials belonging to the actual condition 
    if ismember(actCond.blocktype,'pre-test')
        cfg.trials = ismember(cond,actCond.condition{:}) & ...
            badTrials == 0 & ~catchTrials;
    else
        cfg.trials = ismember(cond,actCond.condition{:}) & ...
            badTrials == 0 & ismember(adaptdir,actCond.adaptdir) & ...
            ~catchTrials;
    end
    ftDataAvg = ft_timelockanalysis(cfg,ftDataDsampl);
    % Getting rid of unnecesary previous cfgs
    ftDataAvg.cfg.previous = [];
    
    % Generating the condition tag
    blockType = actCond.blocktype{:};
    if actCond.aloc < 0
        locStr = 'leftSide';
    else
        locStr = 'rightSide';
    end
    adaptdirStr = actCond.adaptdir{:};
    if isempty(adaptdirStr)
        condStr = [subID,'_',blockType,'_',locStr];
    else
        condStr = [subID,'_',blockType,'_',adaptdirStr,'_',locStr];
    end
    
    %% Saving data
    savePath = fullfile(destDir,['fteeg_ERP_',condStr,'.mat']);
    save(savePath,'ftDataAvg','-v7.3');
    
    ftDataAvg = [];
    
    parfor_progress;
end

% Finalizing progress monitor.
parfor_progress(0);

end


function behavData = loadBehavData(behavSourceDir,fileSpec)

saveDf = cd(behavSourceDir);
fileList = cellstr(ls('*.mat'));
runsFound = regexp(fileList,'.*_run([0-9]+)_.*','tokens','once');
runsFound = [runsFound{:}]';
fileList = fileList(ismember(cellfun(@str2num,runsFound),fileSpec.runid));

for i = 1:numel(fileList)
    behavData(i) = load(fileList{i}); %#ok
end

cd(saveDf);

end