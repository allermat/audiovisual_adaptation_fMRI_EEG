function badTrials = rejecttrials(eegData,behavData,artfData,trigDef,winOfInt,behavSourceDir)
% Checks if trials meet various criteria
% 
% DETAILS:
%   Checks each trial found in the eegData whether: 
%       (1) the response is missing
%       (2) there is a false alarm
%       (3) the response was made by wrong hand (not applicable)
%       (4) there is an EEG artefact
%       (5) there is eyetracker data present for the trial
%       (6) there is an eyeblink
%       (7) there is a saccade
%       (8) the fixation is at the correct location
%       
% INPUT:
%   eegData: epoched fieldtrip data file
%   behavData: sturcture array of behavioural data (each array
%       element contains data from one run)
%   artfData: structure array of detected EEG artefacts (corresponding to
%       the fieldtrip data file)
%   trigDef: table/dataset array of trigger code definitions
%   winOfInt: 1 x 2 vector with the start and end time of time window of 
%       interest with respect to the target event onset
%   behavSourceDir: path to the behavioural source files' folder
% 
% OUTPUT: 
%   badTrials: vector of number of trials x 1 with the following codes
%       0 - good
%       1 - no response
%       2 - false alarm
%       3 - wrong hand (not used)
%       4 - EEG artefact
%       5 - missing eyetracker data
%       6 - eyeblink
%       7 - saccade
%       8 - wrong fixation location
% 

% Parsing input
p = inputParser;

% Input checking functions
checkFtData = @(x) ft_datatype(x,'raw');
checkBehavData = @(x) isstruct(x) && ...
    all(isfield(x,{'run','data'}));
checkArtfData = @(x) isstruct(x) && all(ismember(fieldnames(x),{'artefact_muscle','artefact_visual'}));
checkTrigDef = @(x) isa(x,'dataset') && ...
    all(ismember({'aloc','blocktype','condition'},x.Properties.VarNames));

% Defining input
addRequired(p,'eegData',checkFtData);
addRequired(p,'behavData',checkBehavData);
addRequired(p,'artfData',checkArtfData);
addRequired(p,'trigDef',checkTrigDef);
addRequired(p,'winOfInt',@(x)validateattributes(x,{'numeric'},{'size',[1,2],'finite','nonnan'}));
addRequired(p,'behavSourceDir',@(x) exist(x,'dir'));

% Parsing inputs
parse(p,eegData,behavData,artfData,trigDef,winOfInt,behavSourceDir);

% Assigning input to variables
eegData = p.Results.eegData;
behavData = p.Results.behavData;
artfData = p.Results.artfData;
trigDef = p.Results.trigDef;
winOfInt = p.Results.winOfInt;
behavSourceDir = p.Results.behavSourceDir;

% eeg data properties
nTrials = size(eegData.trial,2);
trialInfo = eegData.trialinfo;
actRun = trialInfo(1,2);
actBlock = trialInfo(1,3);

% Behavioural data
temp = [behavData.run];
runIDs = [temp.id];
behavDataActRun = behavData(runIDs == actRun);
behavDataActBlock = behavDataActRun.data(...
    behavDataActRun.data.block == actBlock,:);
iTiralInBlock = 1;

% EyeLink data
eyeDataActRun = analyse_eyedata(fullfile(behavSourceDir,behavDataActRun.fname),'plot',0);
eyeDataActBlock = eyeDataActRun(behavDataActRun.data.block == actBlock);

% Array for collecting trial rejection info
badTrials = zeros(nTrials,1);

for iTrialInFile = 1:nTrials
    
    % Updating the actual  run, block, and related
    % variables if necessary
    if actRun ~= trialInfo(iTrialInFile,2)
        actRun = trialInfo(iTrialInFile,2);
        % Loading the behavioural data file corresponding to the
        % actual run
        behavDataActRun = behavData(runIDs == actRun);
        % Getting the eyelink data corresponding to the actual run
        eyeDataActRun = analyse_eyedata(fullfile(behavSourceDir, ...
                                                 behavDataActRun.fname),'plot',0);
        % Setting actBlock to NaN so the next conditional will be
        % executed after the run has been changed no matter what
        % the block number was
        actBlock = NaN;
    end
    if actBlock ~= trialInfo(iTrialInFile,3)
        actBlock = trialInfo(iTrialInFile,3);
        behavDataActBlock = behavDataActRun.data(...
            behavDataActRun.data.block == actBlock,:);
        eyeDataActBlock = eyeDataActRun(behavDataActRun.data.block == actBlock);
        iTiralInBlock = 1;
    end

    % Checking behavioural data
    % Sanity check 
    actCond = trialInfo(iTrialInFile,end);
    if behavDataActBlock.aloc(iTiralInBlock) ~= trigDef.aloc(trigDef.condition == actCond) || ...
            ~strcmp(behavDataActBlock.blocktype{iTiralInBlock},trigDef.blocktype{trigDef.condition == actCond})
        error('rejecttrials:trialMismatch', ...
              ['The behavioual data and the EEG trigger' ...
               'does not match']);
    end
    % Missing response
    if isnan(behavDataActBlock.resp(iTiralInBlock)) && ...
            behavDataActBlock.catch_trial(iTiralInBlock)
        badTrials(iTrialInFile) = 1;
        iTiralInBlock = iTiralInBlock + 1;
        continue;
    end
    % False alarm
    if ~isnan(behavDataActBlock.resp(iTiralInBlock)) && ...
            ~behavDataActBlock.catch_trial(iTiralInBlock)
        badTrials(iTrialInFile) = 2;
        iTiralInBlock = iTiralInBlock + 1;
        continue;
    end
    
    % Checking EEG artefacts
    % Is there an EEG artefact?
    if checkartefacts(eegData.sampleinfo(iTrialInFile,:),eegData.time{iTrialInFile},eegData.fsample,winOfInt,artfData)
        badTrials(iTrialInFile) = 4;
        iTiralInBlock = iTiralInBlock + 1;
        continue;
    end
    
    % Saving eyedata
    if eyeDataActBlock(iTiralInBlock) ~= 0
        badTrials(iTrialInFile) = eyeDataActBlock(iTiralInBlock);
    end
    
    iTiralInBlock = iTiralInBlock + 1;
    
end



end


function foundArtefact = checkartefacts(sampleInfo,time,Fs,winOfInt,artfData)
% Checks whether there is any artefact within the trial
% 
% INPUT: 
%   sampleInfo: 1 x 2 vector of the start and end samples of the trial of
%       interest
%   time: 1 x N vector of time values (in seconds) corresponding to the 
%       samples of the trial of interest, where N is the number of samples.
%   Fs: sampling frequency
%   winOfInt: 1 x 2 vector with the start and end time of time window of 
%       interest with respect to the target event onset
%   artfData: structure array of artefact definitions. Each field of the
%       array is an N x 2 matrix of artefact start and end samples, where N
%       is the number of artefacts. Different fields contain different
%       artefact types. 
% 
% OUTPUT:
%   foundArtefact: true if there is an artefact within the window of
%       interest of the trial
% 

foundArtefact = false;

targEventOnsetSample = sampleInfo(1)-(time(1)*Fs);
startSampleWinOfInt = targEventOnsetSample+(winOfInt(1)*Fs);
endSampleWinOfInt = targEventOnsetSample+(winOfInt(2)*Fs);

artfTypes = fieldnames(artfData);

for i = 1:size(artfTypes,1)
    
    actArtfData = artfData.(artfTypes{i});
    
    % Skip if there are no artefacts
    if isempty(actArtfData)
        break;
    end
    
    % Find artefacts which:
    % begin within the window of interest
    crit1 = actArtfData(:,1) >= startSampleWinOfInt & actArtfData(:,1) <= endSampleWinOfInt;
    % end within the winow of interest
    crit2 = actArtfData(:,2) >= startSampleWinOfInt & actArtfData(:,2) <= endSampleWinOfInt;
    % cover the whole window of interest
    crit3 = actArtfData(:,1) <= startSampleWinOfInt & actArtfData(:,2) >= endSampleWinOfInt;
    
    % If there is an artefact which fall under either criteria, mark it and
    % break. 
    foundArtefact = any(crit1 | crit2 | crit3);
    if foundArtefact
        break;
    end
    
end

end

