function preproc_COMM(subID,varargin)
% AVRECAL project EEG preprocessing common stage

%% Parsing input, checking matlab
p = inputParser;

addRequired(p,'subID',@ischar);
addParameter(p,'overwrite',false,@(x) validateattributes(x, ...
                                                  {'logical'},{'scalar'}));
parse(p,subID,varargin{:});

subID = p.Results.subID;
overwrite = p.Results.overwrite;

% Loading files specifying parameters
trigDef = load(fullfile(get_path('project'),'experiment','trigger_LUT.mat'));
trigDef = trigDef.trigger_LUT;
setupSpec = load(fullfile(get_path('project'),'analysis','EEG','setup_spec.mat'));
setupSpec = setupSpec.setup_spec;
subjectSpec = subject_info;
subjectSpec = subjectSpec(ismember({subjectSpec.id}',subID));
if isempty(subjectSpec)
    error('preproc_COMM:invalidInput', ...
          'Invalid subject ID');
end
if isempty(fieldnames(subjectSpec.EEG))
    error('preproc_COMM:missingData', ...
          'Subject specification does not contain EEG sesssions');
end

analysisDir = fullfile(get_path('project'),'data',subID,'EEG','preproc_data','COMM');
if ~exist(analysisDir,'dir');
    mkdir(analysisDir);
end

tag = 'COMM_';

for iSession = 1:numel(subjectSpec.EEG.session)
    
    actSessionFolderName = subjectSpec.session_folder_name{subjectSpec.EEG.session(iSession)};
    dataDir = fullfile(get_path('project'),'data',subID,'EEG','raw_data',actSessionFolderName);
    if ~exist(dataDir,'dir')
        error('The specified data folder does not exist!');
    end
    
    sourceFileNamesActSession = {subjectSpec.EEG.filespec(iSession).file.name}';
    
    for iFileActSession = 1:size(sourceFileNamesActSession,1)
        
        %% Checking if the processed data are already present
        artfResultFileName = ['artf_',tag,sourceFileNamesActSession{iFileActSession},'.mat'];
        eegResultFileName = ['fteeg_',tag,sourceFileNamesActSession{iFileActSession},'.mat'];
        if exist(fullfile(analysisDir,artfResultFileName),'file') && ...
                exist(fullfile(analysisDir,eegResultFileName),'file')&&...
                ~overwrite;
            warning('Skipping %s as it has been already processed! ',sourceFileNamesActSession{iFileActSession});
            continue;
        end
        
        %% Reading raw data
        cfg = struct();
        cfg.datafile = fullfile(dataDir,[sourceFileNamesActSession{iFileActSession},'.eeg']);
        cfg.headerfile = fullfile(dataDir,[sourceFileNamesActSession{iFileActSession},'.vhdr']);
        cfg.channel = 'all';
        ftDataRaw = ft_preprocessing(cfg);
        
        %% High-pass filtering
        cfg = struct();
        cfg.hpfilter = 'yes';
        cfg.hpfreq = 0.1;
        cfg.hpfiltord = 4;
        cfg.hpfiltdir = 'twopass';
        cfg.plotfiltresp = 'yes';
        fprintf('\nFILTERING: Highpass filter of order %d, half power frequency %d\n\n',cfg.hpfiltord,cfg.hpfreq);
        ftDataHp = ft_preprocessing(cfg,ftDataRaw);
        
        %% Artefact detection and visual check

        % Coarse visual pass on the data just to check if any channels
        % should be rejected
        cfg = struct();
        cfg.headerfile = fullfile(dataDir,[sourceFileNamesActSession{iFileActSession},'.vhdr']);
        cfg.trialdef.eventtype = 'Stimulus';
        cfg.trialdef.faketrllength = 30;
        cfg.trialdef.prestim = 0.1;
        cfg.trialdef.poststim = 1;
        cfg.trialdef.trigdef = trigDef;
        cfg.trialdef.fileSpec = subjectSpec.EEG.filespec(iSession).file(iFileActSession);
        cfg.trialdef.trig_visonset_corr = setupSpec.trig_visonset_corr_eeg;
        cfg.trialfun = 'ft_trialfun_artefactdetection';
        cfg = ft_definetrial(cfg);
        trlArtf = cfg.trl;
        
        % Epoching data for channel check
        cfg = struct();
        cfg.trl = trlArtf;
        ftDataEpArtf = ft_redefinetrial(cfg,ftDataHp);
        
        % Inspect data
        cfg = struct();
        cfg.viewmode  = 'vertical';
        cfg.continuous = 'no';
        cfg.channel = 'all';
        cfg = ft_databrowser(cfg,ftDataEpArtf); %#ok
        
        % Marking bad channels
        channelsList = ftDataHp.label;
        selection = listdlg('ListString',channelsList, ...
                            'SelectionMode', 'multiple', 'ListSize', [200 400],...
                            'PromptString','Choose bad channels');
        % Interpolate bad channels if necessary
        if ~isempty(selection)
            cfg = struct();
            cfg.method = 'spline';
            cfg.badchannel = channelsList(selection);
            cfg.layout = 'eeg1010.lay';
            ftDataHp = ft_channelrepair(cfg,ftDataHp);
            % making sure that the channel labels are in the
            % same order as before interpolation
            [~,idx] = ismember(channelsList,ftDataHp.label);
            ftDataHp.label = ftDataHp.label(idx,:);
            ftDataHp.trial = cellfun(@(x) x(idx,:), ...
                                     ftDataHp.trial,'UniformOutput',false);
        end
        
        % Re-referencing to average reference
        cfg = struct();
        cfg.reref = 'yes';
        cfg.refchannel = 'all';
        ftDataReref = ft_preprocessing(cfg,ftDataHp);
        % Clearing unnecessary previous dataset
        ftDataHp = []; %#ok<NASGU>
        
        % Automatic artefact detection
        % Marking trials for artefact detection, now with a
        % finer grain epoching
        cfg = struct();
        cfg.headerfile = fullfile(dataDir,[sourceFileNamesActSession{iFileActSession},'.vhdr']);
        cfg.trialdef.eventtype = 'Stimulus';
        cfg.trialdef.faketrllength = 15;
        cfg.trialdef.prestim = 0.1;
        cfg.trialdef.poststim = 1;
        cfg.trialdef.trigdef = trigDef;
        cfg.trialdef.fileSpec = subjectSpec.EEG.filespec(iSession).file(iFileActSession);
        cfg.trialdef.trig_visonset_corr = setupSpec.trig_visonset_corr_eeg;
        cfg.trialfun = 'ft_trialfun_artefactdetection';
        cfg = ft_definetrial(cfg);
        trlArtf = cfg.trl;
        
        % Epoching data for automatic artefact detection
        cfg = struct();
        cfg.trl = trlArtf;
        ftDataEpArtf = ft_redefinetrial(cfg,ftDataReref);
        
        % - MUSCLE ARTEFACTS -
        cfg = struct();
        cfg.continuous = 'no';
        % channel selection, cutoff and padding
        cfg.artfctdef.zvalue.channel = 'all'; %{'FT7','FT8','AF7','AF8','TP7','TP8','PO7','PO8'};
        cfg.artfctdef.zvalue.cutoff = subjectSpec.EEG.filespec(iSession).file(iFileActSession).cutoff_zval;
        cfg.artfctdef.zvalue.trlpadding = 0;
        cfg.artfctdef.zvalue.fltpadding = 0;
        cfg.artfctdef.zvalue.artpadding = 0.1;
        % algorithmic parameters
        cfg.artfctdef.zvalue.bpfilter = 'yes';
        cfg.artfctdef.zvalue.bpfreq = [110 140];
        cfg.artfctdef.zvalue.bpfilttype = 'but';
        cfg.artfctdef.zvalue.bpfiltord = 9;
        cfg.artfctdef.zvalue.hilbert = 'yes';
        cfg.artfctdef.zvalue.boxcar = 0.5;
        % make the process interactive
        cfg.artfctdef.zvalue.interactive = 'no';
        
        [~,artefact_muscle] = ft_artifact_zvalue(cfg,ftDataEpArtf);
        
        % - JUMP ARTEFACTS -
        %         cfg = struct();
        %         % required fields
        %         cfg.continuous                  = 'no';
        %         cfg.artfctdef.zvalue.channel    = 'all';
        %         cfg.artfctdef.zvalue.cutoff     = 60; % µV (let's assume that if there are jumps, they are pretty big)
        %         cfg.artfctdef.zvalue.trlpadding = 0;  % add a bit of data on both sides
        %         cfg.artfctdef.zvalue.fltpadding = 0;  % for trial and filter padding
        %         cfg.artfctdef.zvalue.artpadding = 0;  % delete time period on both sides
        %         % algorithmic parameters
        %         cfg.artfctdef.zvalue.cumulative    = 'yes';
        %         cfg.artfctdef.zvalue.medianfilter  = 'yes'; % preserves jumps
        %         cfg.artfctdef.zvalue.medianfiltord = 9;
        %         cfg.artfctdef.zvalue.absdiff       = 'yes';
        %         cfg.artfctdef.zvalue.interactive   = 'no';
        %
        %         [~,artefact_jump] = ft_artifact_zvalue(cfg,dataEpArtf);
        
        %% Inspect detected artifacts
        cfg = struct();
        cfg.viewmode  = 'vertical';
        cfg.continuous = 'no';
        cfg.channel = 'all';
        cfg.artfctdef.muscle.artifact = artefact_muscle;
        %         cfg.artfctdef.jump.artifact = artefact_jump;
        cfg = ft_databrowser(cfg,ftDataEpArtf);
        
        %% Saving artefact data
        fprintf('\n\nSaving artefact data...\n\n');
        % artefacts
        artefact_muscle = cfg.artfctdef.muscle.artifact; %#ok<NASGU>
        artefact_visual = cfg.artfctdef.visual.artifact; %#ok<NASGU>
        
        save(fullfile(analysisDir,artfResultFileName),'artefact_muscle','artefact_visual','-v7.3');

        
        %% Saving EEG data
        fprintf('\n\nSaving EEG data...\n\n');
        save(fullfile(analysisDir,eegResultFileName),'ftDataReref','-v7.3');
        
    end
end

end
