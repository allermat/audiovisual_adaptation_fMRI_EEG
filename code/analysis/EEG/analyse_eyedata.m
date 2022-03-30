function [badTrials, values] = analyse_eyedata(behavfname, varargin)
% Quick analysis for assessing eyelink data quality of a given run

if ~mod(nargin, 2)
    error('variable input arguments should come in pairs')
end
S = struct('fixation', 'true', 'plot', 1, 'block', [], 'rawread', 0, 'esp', 0);
for i=1:2:numel(varargin)
    S.(varargin{i}) = varargin{i+1};
end
if ~ismember(S.fixation, {'relative' 'true'})
   error('fixation should be true or relative'); % cluster 
end

if ~exist('behavfname', 'var')
    error('file name with full path should be given');
end

% Load behavioural data file
load(behavfname, 'data', 'setup', 'subj');
behavData = data;
[pathstr, elfname, ext] = fileparts(setup.eyelink.edfFileName);

%% Preprocessing data
load(fullfile(get_path('project'), 'experiment', 'trigger_LUT.mat'));

% Preprocessing eyetracking data
% subj.session = 1;
subjInfo = subject_info;
subjInfo = subjInfo(ismember({subjInfo.id},subj.id));
if isempty(subjInfo)
    error('analyse_eyedata:invalidSubject','Subject not found');
end
if isfield(subjInfo, 'session_folder_name') && ~isempty(subjInfo.session_folder_name)
    subj.session_folder_name = subjInfo.session_folder_name;
end
eyefolder = fullfile(get_path('project'), get_folder(subj, 'r', 'eyetracking'));
edfname = [elfname '_eyedata.mat'];
if ~S.rawread && exist(fullfile(eyefolder, edfname), 'file')
    load(fullfile(eyefolder, edfname));
else
    ascFilePath = fullfile(eyefolder, [elfname '.asc']);
    if exist(ascFilePath, 'file')
        eyeData = preproceye(ascFilePath, trigger_LUT, setup);
        save(fullfile(eyefolder, edfname), 'eyeData');
    else
        fprintf('%s\n', ascFilePath);
        error('asc file not found, probably due to conversion error.');
    end
end

% Filter to a specific block
if ~isempty(S.block)
    blockStartId = find(ismember({eyeData.event.type},'Stimulus') & ...
        ismember({eyeData.event.value},trigger_LUT.trig_eye(ismember(trigger_LUT.type,'blockstart'))));
    blockEndId = find(ismember({eyeData.event.type},'Stimulus') & ...
        ismember({eyeData.event.value},trigger_LUT.trig_eye(ismember(trigger_LUT.type,'blockend'))));
    
    if numel(blockStartId) < S.block
        warning('requested block is not present in data');
        [badTrials, values] = deal(NaN);
        return;
    end
    eyeData.event = eyeData.event(blockStartId(S.block):blockEndId(S.block));
end

% Start and end sample of events
evStartSamples = [eyeData.event.sample]';
evEndSamples = evStartSamples+[eyeData.event.duration]';
targEvStartSamples = evStartSamples(ismember({eyeData.event.type},'Stimulus') & ...
    ismember({eyeData.event.value},trigger_LUT.trig_eye(ismember(trigger_LUT.type,'stim'))));


%% Saccade, RT and fixation criteria

% Load criteria from file
load(fullfile(get_path('project'), 'experiment', 'behav_criteria.mat'), 'eyetracking');

% Saccade criteria
criteria.saccade = eyetracking.saccade;

% Reaction time criteria
criteria.RT = eyetracking.RT;

% Fixation criteria
fixcoord = cellfun(@str2num, {eyeData.event(ismember({eyeData.event.type}, 'Fixation')).value}', 'un', 0);
fixcoord = cat(1, fixcoord{:});
medianX = median(fixcoord(:,1)); % median x coordianate of fixation in data
medianY = median(fixcoord(:,2)); % median y coordianate of fixation in data
[trueX, trueY] = RectCenter(setup.eyelink.rect); % true x,y coordinates of fixation
if ~isPointInCircle([trueX trueY], [medianX, medianY deg2pix(eyetracking.fixation.offset, setup)])
    offset = sqrt((medianX - trueX)^2 + (medianY - trueY)^2);
    warning('True fixation point is off with %.2f degrees', pix2deg(offset, setup));
end
criteria.fixation.r = deg2pix(eyetracking.fixation.radius, setup); % allowed radius around the fixation cross
criteria.fixation.x = trueX;
criteria.fixation.y = trueY;

% Check saccade data
value = {eyeData.event(ismember({eyeData.event.type},'Saccade')).value}';
value = cell2mat(cellfun(@str2num, value, 'un', 0));
if S.esp && size(value, 2) ~= 4
    warning('saccade do not contain end of position values');
end

%% Plot fixation points

if S.plot
    fig = plotfixation(eyeData, behavData.RT, criteria, setup, targEvStartSamples, evStartSamples, evEndSamples);
    saveas(fig, fullfile(eyefolder, [elfname '_fixation']),'png');
end

if strcmp(S.fixation, 'relative')
    criteria.fixation.x = medianX;
    criteria.fixation.y = medianY;
end

%% Marking and plotting bad trials

% Last column of trialinfo gives the following info:
% 0 - good
% 5 - missing eyetracker data
% 6 - eyeblink
% 7 - saccade
% 8 - wrong fixation location
[badTrials, values] = checktrials(eyeData, behavData, criteria, targEvStartSamples, evStartSamples, evEndSamples);

% Plot
if S.plot
    fig = plotbadtrials(badTrials);
    saveas(fig, fullfile(eyefolder, [elfname '_quickanal']),'png');
end
    
end

%% --------------- Preprocessing functions ------------------

function dataEye = preproceye(ascFilePath,trigDef,setupSpec)

% stimTriggers = trigDef.trig_eye(trigDef.type == 'stim');
% trig_visonset_corr_eyelink = setupSpec.trig_visonset_corr_eyelink;
stimTriggers = trigDef.trig_eye(ismember(trigDef.type, 'stim'));
trig_visonset_corr_eyelink = setupSpec.eyelink.delay;

dataEye = struct('event',[],'Fs',[]);

[event,hdr] = read_eyelink_event(ascFilePath);

evValues = {event.value}';
evTypes = {event.type}';
evStartSamples = [event.sample]';
% Finding the onset samples of stimulus triggers
targEvStartSamples = evStartSamples(ismember(evTypes,'Stimulus') & ...
    ismember(evValues,stimTriggers));
% Correcting stimulus triggers for trigger-visual onset asynchrony
targEvStartSamples = targEvStartSamples+(trig_visonset_corr_eyelink*hdr.Fs);
% Saving corrected values into the original structure
evStartSamples(ismember(evTypes,'Stimulus') & ...
    ismember(evValues,stimTriggers)) = targEvStartSamples;
evStartSamples = num2cell(evStartSamples);
[event.sample] = evStartSamples{:};

dataEye.event = event;
dataEye.Fs = hdr.Fs;

end

function idEvent = selecteyeevents(actTrial, RT, winOfInt, targEvStartSamples, evStartSamples, evEndSamples, Fs)
% Selects eyetracker events corresponding to a particular trial
% 
% INPUT: 
%   actTrial           = serial number of the trial of interest within its session
%   RT                 = behavioural reaction time data
%   winOfInt           = 1 x 2 vector with the start and end time of time window of 
%                           interest with respect to the target event onset
%   targEvStartSamples = starting sample of target stimulus events
%   evStartSamples     = starting sample of all events
%   evEndSamples       = ending sample of all events
%   Fs                 = eyedata sampling frequency 
% 
% OUTPUT:
%   idEvent = index of events corresponding to the trial of interest.
%               Empty array if no events were found.
% 

% If the trial is not found return an empty id
if ~ismember(1:size(targEvStartSamples,1),actTrial)
    idEvent = [];
    return;
end

% Find window of interest around target event
startSampleWinOfInt = targEvStartSamples(actTrial)+(winOfInt(1)*Fs);
if numel(winOfInt) == 2
    endSampleWinOfInt = targEvStartSamples(actTrial)+(winOfInt(2)*Fs);
else % window up to response if no end of window
    endSampleWinOfInt = targEvStartSamples(actTrial)+(RT(actTrial)*Fs);
end

% Find events which:
% begin within the window of interest
idEvent(:,1) = evStartSamples >= startSampleWinOfInt & evStartSamples <= endSampleWinOfInt;
% end within the winow of interest
idEvent(:,2) = evEndSamples >= startSampleWinOfInt & evEndSamples <= endSampleWinOfInt;
% cover the whole window of interest
idEvent(:,3) = evStartSamples <= startSampleWinOfInt & evEndSamples >= endSampleWinOfInt;

end

%% -------------- Criteria checking functions ------------------

function [badTrials, values] = checktrials(eyeData, behavData, criteria, targEvStartSamples, evStartSamples, evEndSamples)

nTrials = size(behavData,1);
% Array for collecting rejection info
[badTrials, values] = deal(zeros(nTrials,1));

for actTrialInSession = 1:nTrials
    % Is there a response? 
%     if isnan(behavData.resp(actTrialInSession)) && behavData.catch_trial(actTrialInSession) % behavData.iTrialInSession == 
%         badTrials(actTrialInSession) = 1;
%         continue;
%     end
    
    % Is the response too early?
%     if behavData.RT(actTrialInSession) < criteria.RT.min && any(behavData.catch_trial) % behavData.iTrialInSession == respTime
%         badTrials(actTrialInSession) = 2;
%         continue;
%     end
    
    % Choose events around the window of interest of the stimulus
    if ~isempty(strfind(behavData.blocktype{actTrialInSession}, 'test')) && ~behavData.catch_trial(actTrialInSession)
        idEvent = selecteyeevents(actTrialInSession, behavData.RT, criteria.RT.window.long, targEvStartSamples, evStartSamples, evEndSamples, eyeData.Fs);
    else
        idEvent = selecteyeevents(actTrialInSession, behavData.RT, criteria.RT.window.short, targEvStartSamples, evStartSamples, evEndSamples, eyeData.Fs);
    end
    eyeEventsActTrial = eyeData.event(any(idEvent, 2));
    if isempty(eyeEventsActTrial)
%         warning('No eyetracking data for trial %d',actTrialInSession);
        badTrials(actTrialInSession) = 5;
        values(actTrialInSession) = -1;
        continue;
    end
    
    % Is there a blink? 
    if any(ismember({eyeEventsActTrial.type},'Blink'))
        badTrials(actTrialInSession) = 6;
        continue;
    end
    
    % Is there a saccade? 
    if any(ismember({eyeEventsActTrial.type},'Saccade'))
        
        saccades = eyeEventsActTrial(ismember({eyeEventsActTrial.type},'Saccade'));
        % Are any of the saccades above threshold?
        if checksaccades(saccades, criteria.saccade, eyeData.Fs)
            badTrials(actTrialInSession) = 7;
            saccendcoord = cell2mat(cellfun(@str2num, {saccades.value}', 'un', 0));
            values(actTrialInSession) = criteria.fixation.x - mean(saccendcoord(:,1));
            continue;
        end
        
    end
    
    % Is fixation average location correct?
    if all(~ismember({eyeEventsActTrial.type},'Fixation'))
        badTrials(actTrialInSession) = 8;
        continue;
    else
        fixations = eyeEventsActTrial(ismember({eyeEventsActTrial.type},'Fixation'));
        % Are any of the saccades above threshold?
        if checkfixations(fixations, criteria.fixation)
            badTrials(actTrialInSession) = 8;
        end
        fixcoord = cell2mat(cellfun(@str2num, {fixations.value}', 'un', 0));
        values(actTrialInSession) = criteria.fixation.x - mean(fixcoord(:,1));
    end
    
end

end


function aboveThreshold = checksaccades(saccades, saccadeCriteria, Fs)
% Checks whether saccade parameters are above the given threshold
% 
% INPUT: 
%   saccades: structure array containing saccade events
%   thresholds: 1 x 3 vector of the threshold values for defining a saccade
%               thresholds(1) = minimum amplitude in degrees
%               thresholds(2) = minimum peak velocity in degrees/seconds
%               thresholds(3) = minimum duration in seconds
%   Fs: sampling frequency
% 
% OUTPUT:
%   aboveThreshold: true if any of the saccade events are above all three 
%       thresholds. 
% 

if any(~ismember({saccades.type},'Saccade'))
    error('Event of wrong type passed to checksaccades!');
end

temp = regexp({saccades.value},'([0-9.e+]+)\s+([0-9.e+]+)','tokens');
if isempty(temp)
    error('Wrong saccade value format!');
end
temp = [temp{:}]';
temp = vertcat(temp{:});
temp = cell2mat(cellfun(@str2num,temp,'UniformOutput',false));
ampl = temp(:,1);
pv = temp(:,2);
% Converting duration from samples to seconds
dur = [saccades.duration]'./Fs;

if ampl >= saccadeCriteria.ampl % & pv >= saccadeCriteria.pv & dur >= saccadeCriteria.dur)
    aboveThreshold = true;
else
    aboveThreshold = false;
end

end

function failedCriteria = checkfixations(fixations, fixCriteria)
% Checks whether fixation parameters meet the given criteria
% 
% INPUT: 
%   fixations: structure array containing fixation events
%   criteria: vector of the criteria for accepting a fixation event
%               criteria(1) = fixation cross x coordinate (in pixels)
%               criteria(2) = fixation cross y coordinate (in pixels)
%               criteria(3) = radius of the accepted area around the
%                             fixation cross (in pixels)
% 
% OUTPUT:
%   failedCriteria: true if any of the fixation events fail to meet the
%       criteria
% 

if any(~ismember({fixations.type},'Fixation'))
    error('Event of wrong type passed to checkfixations!');
end

temp = regexp({fixations.value},'([0-9.e+]+)\s+([0-9.e+]+)','tokens');
if isempty(temp)
    error('Wrong saccade value format!');
end
temp = [temp{:}]';
temp = vertcat(temp{:});
if isempty(temp)
    failedCriteria = true;
    return;
else
    temp = cell2mat(cellfun(@str2num,temp,'UniformOutput',false));
end
% This should produce an N x 2 matrix, where N is the number of fixation
% events and the two columns are the average x and y positions 
% respectively. 

if any(~isPointInCircle(temp,repmat([fixCriteria.x fixCriteria.y fixCriteria.r], size(temp,1),1)))
    failedCriteria = true;
else
    failedCriteria = false;
end

end


%% ------------------ Plotting functions -----------------

function [fig, C] = plotfixation(eyeData, RT, criteria, setup, targEvStartSamples, evStartSamples, evEndSamples)

nTrials = size(RT, 1);

fig = figure;
hold on;

fixcoord = [];

for actTrialInSession = 1:nTrials
    % Choose events around the window of interest of the stimulus
    idEvent = selecteyeevents(actTrialInSession, RT, criteria.RT.window.short, targEvStartSamples, evStartSamples, evEndSamples, eyeData.Fs);
    eyeEventsActTrial = eyeData.event(any(idEvent, 2));

    if ~isempty(eyeEventsActTrial)
        idfixation = ismember({eyeEventsActTrial.type},'Fixation');
        
        if any(idfixation)
            fixations = eyeEventsActTrial(idfixation);
            temp = regexp({fixations.value},'([0-9.e+]+)\s+([0-9.e+]+)','tokens');
            temp = [temp{:}]';
            temp = vertcat(temp{:});
            if ~isempty(temp)
                fixcoord = [fixcoord; cell2mat(cellfun(@str2num,temp,'UniformOutput',false))];
            end
        end
    end
end

% C = find_fixation_clusters(fixcoord, fixcriteria);
if isempty(fixcoord)
    warning('no fixation found in data')
else
    plot(fixcoord(:,1), fixcoord(:,2), 'bo');
end
% find_fixation_clusters(fixcoord, fixcriteria, C);

xlim([setup.eyelink.rect(1) setup.eyelink.rect(3)]);
ylim([setup.eyelink.rect(2) setup.eyelink.rect(4)]);

fixCriteria = criteria.fixation;
plot([-fixCriteria.r fixCriteria.r]+fixCriteria.x, [fixCriteria.y fixCriteria.y], 'c');
plot([fixCriteria.x fixCriteria.x], [-fixCriteria.r fixCriteria.r]+fixCriteria.y, 'c');

end

function fig = plotbadtrials(badTrials)

id = [1 6:9];
exclLevelsAll = 0:8;
exclLevelsBad = [NaN,1:8];

freqOfLevelsAll = calcfreqoflevels(badTrials,exclLevelsAll(id));
freqOfLevelsBad = calcfreqoflevels(badTrials(badTrials ~= 0),exclLevelsBad(id));

fig = figure();
set(gcf,'Units','normalized','Position',[0.7,0.3,0.225,0.3],...
    'PaperPositionMode','auto');
    
h = bar(1:2,[freqOfLevelsAll;freqOfLevelsBad],'stacked');
ylabel('Proportion of trials');
ylim([0,1]);
set(gca,'XTick',[],'Box','off');
axesPosition = get(gca,'Position');
axesPosition(3) = 0.5*axesPosition(3);
set(gca,'Position',axesPosition);
if freqOfLevelsAll(1) < 1
    yMaxRightYaxis = 1-freqOfLevelsAll(1);
else
    yMaxRightYaxis = 1;
end
hNewAxes = axes('Position',axesPosition,'Color','none','YLim',[0 yMaxRightYaxis],...
                'YAxisLocation','right','XTick',[],'Box','off');
legendstr = {'good','no resp','early resp','wrong hand','EEG artf','no EYE data','blink','saccade','wrong fixation'};
hLeg = legend(h,legendstr(id),'Location','SouthEastOutside');
% Offsetting legend position
legendPosition = get(hLeg,'Position');
legendPosition(1) = axesPosition(1)+(axesPosition(3)*1.2);
set(hLeg,'Position',legendPosition);

end


function out = calcfreqoflevels(X,levels)

resp = X(:,1);
nLevels = numel(levels);
out = zeros(1,nLevels);

for i = 1:nLevels
    out(i) = sum(resp == levels(i))/size(resp,1);
end

end


function C = find_fixation_clusters(fixcoord, fixcriteria, C)

if exist('C', 'var')
    while 1
        str = input('Enter the id of cluster centroids, then quit (q): ', 's');
        if strcmp(str, 'q')
            break;
        else
            k = str2num(str);
            viscircles(C(k,:), fixcriteria.r, 'EdgeColor', 'r');
        end
    end
else
    while 1
        plot(fixcoord(:,1), fixcoord(:,2), 'bo');
        if exist('C', 'var')
            viscircles(C, repmat(fixcriteria.r, k, 1), 'EdgeColor', 'r');
        end
        
        str = input('Enter the number of clusters or quit (q): ', 's');
        cla;
        if strcmp(str, 'q')
            break;
        else
            k = str2num(str);
        end
        [IDX, C] = kmeans(fixcoord, k);
        xlabel('x coord');
        ylabel('y coord');
    end
    C
    
    % eucD = pdist(fixcoord,'euclidean');
    % Z = linkage(eucD,'average');
    % color = {'r' 'g' 'm' 'y' 'c'}
    % cla
    % c = cluster(Z, 'maxclust', 11)
    % for kk=1:11
    % viscircles(mean(fixcoord(c==kk,:), 1), fixation.r, 'EdgeColor', 'r');
    % end
    % plot(fixcoord(:,1), fixcoord(:,2), 'bo');
end

end