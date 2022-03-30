function preproc_MVPA(subID,varargin)
%% Parsing input, checking matlab
p = inputParser;

addRequired(p,'subID',@ischar);
addParameter(p,'overwrite',false,@(x) validateattributes(x,{'logical'},{'scalar'}));
parse(p,subID,varargin{:});

subID = p.Results.subID;
overwrite = p.Results.overwrite;

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
destDir = fullfile(get_path('project'),'data',subID,'EEG','preproc_data','MVPA');
if ~exist(destDir,'dir');
    mkdir(destDir);
end

fprintf('\n\nProcessing files...\n');

for iSession = 1:numel(subjectSpec.EEG.session)
    
    actSessionFolderName = subjectSpec.session_folder_name{subjectSpec.EEG.session(iSession)};
    rawDataDir = fullfile(get_path('project'),'data',subID,'EEG','raw_data',actSessionFolderName);
    if ~exist(rawDataDir,'dir')
        error('The specified data folder does not exist!');
    end
    behavSourceDir = fullfile(get_path('project'),'data',subID, ...
                              'EEG','behavioural data',actSessionFolderName);
    
    actSessionStr = num2str(subjectSpec.EEG.session(iSession));
    destFilePath = fullfile(destDir,['fteeg_MVPA_',subID,'_session',actSessionStr,'.mat']);
    if exist(destFilePath,'file') && ~overwrite
        warning('Skipping session %s as it has been already processed! ',actSessionStr);
        continue;
    end
    
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
        
        % if all(~ismember(subjectSpec.EEG.posttest.radapt.runid{iSession},actFileSpec.runid))
        %     warning(['This file does not contain blocks of required ' ...
        %              'type, skipping. ']);
        %     continue;
        % end
        
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
        cfg.trialdef.blocktype = {'pre-test','post-test'};%{'r-adaptation'};
        cfg.trialdef.prestim = 0.1;
        cfg.trialdef.poststim = 0.6;
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
    
    % Cleaning up result array from empty elements
    resultFilesActSession(cellfun(@isempty,resultFilesActSession)) = [];
    
    %% Merging files for the same day if necessary
    if numel(resultFilesActSession) > 1
        ftDataClean = ft_appenddata([],resultFilesActSession{:});
    else
        ftDataClean = resultFilesActSession{1}; %#ok<*NASGU>
    end
    
    resultFilesActSession = [];
    
    %% Saving data
    save(destFilePath,'ftDataClean','-v7.3');
    
    ftDataClean = [];
    
end

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