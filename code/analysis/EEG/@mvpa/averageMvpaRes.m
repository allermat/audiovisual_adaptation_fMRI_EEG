function mvparesObj = averageMvpaRes(I)
% Method for averaging mvpares objects
% 
% USAGE:
%   mvparesObj = averageMvpaRes(I)
% INPUT:
%   I (struct): structure of inputs. Required fields:
%       pathAveragedFile: full path for saving the averaged file
%       pathFilesToAverage: cell array of full paths to the to be averaged
%           files
% OUTPUT:
%   mvparesObj (object): averaged mvpares object

% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

% Parsing input
p = inputParser;
wrapexist = @(x) exist(x,'file');
requiredVars = {'pathAveragedFile','pathFilesToAverage'};
addParameter(p,'pathAveragedFile','');
addParameter(p,'pathFilesToAverage','',@(x) all(cellfun(wrapexist,x)));

parse(p,I);

pathAveragedFile = p.Results.pathAveragedFile;
pathFilesToAverage = p.Results.pathFilesToAverage;

if any(ismember(requiredVars,p.UsingDefaults))
    error('All required parameters must be specified!');
end

% Loading objects from disk
objList = cell(size(pathFilesToAverage));
for i = 1:numel(pathFilesToAverage)
    objList{i} = mvpares(pathFilesToAverage{i});
end
% Checking objects
checkObj = @(x) isa(x,'mvpares') && x.isvalid;
if any(~cellfun(checkObj,objList))
    error('mvpa:averageMvpaRes:invalidInput',...
        'All inputs should be valid instances of mvpares class.');
end
if numel(objList) < 2
    error('mvpa:averageMvpaRes:singleInput',...
        'There is just one input object, can''t average, returning.');
end
avgType = checkInfos(objList);

% Looking for data available for averaging in the objects
dataFields = cellfun(@who,objList,'UniformOutput',false);
if strcmp(avgType,'acrossSub')
    averagedDataNames = {'gen_neuroMetrFun','gen_perfEstimates',...
                        'gen_recalIndex'};
    % averagedDataNames = {'gen_neuroMetrFun','gen_neuroMetrFun_ffx',...
                        % 'gen_perfEstimates','gen_recalIndex', ...
                        % 'gen_recalIndex_ffx'};
else
    averagedDataNames = {'gen_neuroMetrFun','gen_perfEstimates',...
                        'gen_recalIndex'};
end

averages = cell(size(averagedDataNames));

for i = 1:numel(averagedDataNames)
    
    switch averagedDataNames{i}
        case 'gen_neuroMetrFun'
            indivDataToCheck = 'gen_neuroMetrFun';
        case 'gen_neuroMetrFun_ffx'
            indivDataToCheck = 'gen_neuroMetrFun';
        case 'gen_perfEstimates'
            indivDataToCheck = 'gen_perfEstimates';
        case 'gen_recalIndex'
            indivDataToCheck = 'gen_recalIndex';
        case 'gen_recalIndex_ffx'
            indivDataToCheck = 'gen_recalIndex';
    end
    
    if any(cellfun(@ismember,repmat({indivDataToCheck},size(dataFields)),dataFields))
        if any(~cellfun(@ismember,repmat({indivDataToCheck},size(dataFields)),dataFields))
            warning('mvpa:averageMvpaRes:missingData',...
                ['Not all of the input objects have the ''%s'' ',...
                'field. This field will not be averaged.'],averagedDataNames{i});
        else
            averages{i} = averageIndividualData(objList,averagedDataNames{i},avgType);
        end
    end
end

if any(~cellfun(@isempty,averages))
    keepIdx = ~cellfun(@isempty,averages);
    % Removing fields which are not averaged
    averages = averages(keepIdx);
    fieldsAveraged = averagedDataNames(keepIdx);
    avgStruct = cell2struct(averages,fieldsAveraged,2);
    % Adding info field to the averaged dataset
    avgInfo = objList{1}.getInfo;
    infoFieldsToKeep = {'cv_scheme','gen_cond','gen_data_file','gen_data_fs',...
        'gen_label','gen_time','sc_method','subID','svm_type','tr_cond',...
        'tr_data_file','tr_data_fs','tr_label','tr_method','tr_timePoints','timeBinRegFun'};
    infoFields = fieldnames(avgInfo);
    infoFieldsToRemove = infoFields(~ismember(infoFields,infoFieldsToKeep));
    avgInfo = rmfield(avgInfo,infoFieldsToRemove);
    if strcmp(avgType,'acrossSub')
        avgInfo.subID = 'group';
    end
    avgInfo.sourceFiles = pathFilesToAverage;
    avgInfo = orderfields(avgInfo);
    avgStruct.info = avgInfo;
    avgStruct = orderfields(avgStruct); %#ok<NASGU>
    % Saving averaged dataset
    save(pathAveragedFile,'-struct','avgStruct','-v7.3');
    % Loading mvpares object from the saved dataset
    mvparesObj = mvpares(pathAveragedFile);
else
    warning('mvpa:averageMvpaRes:requestedDataNotAvailable',...
        'Could not find any data to average. Retruning ');
    mvparesObj = [];
end

end

function avgType = checkInfos(objList)

infos = cellfun(@getInfo,objList);
% Determining if the average is within or across subjects
if numel(unique({infos.subID})) == size(infos,1)
    avgType = 'acrossSub';
elseif numel(unique({infos.subID})) == 1
    avgType = 'withinSub';
else
    error('mvpa:averageMvpaRes:checkInfos:inconsistentMvpaResults',...
        ['For across subject averaging each subject must have only ' ...
         'one dataset. ']);
end

% Checking for inconsistencies
if numel(unique({infos.cv_scheme})) > 1
    error('mvpa:averageMvpaRes:checkInfos:inconsistentMvpaResults',...
        'The cv_scheme is not consistent across mvpa results for averaging.');
elseif numel(unique({infos.gen_cond})) > 1
    error('mvpa:averageMvpaRes:checkInfos:inconsistentMvpaResults',...
        'The gen_cond is not consistent across mvpa results for averaging.');
elseif numel(unique({infos.gen_label})) > 1
    error('mvpa:averageMvpaRes:checkInfos:inconsistentMvpaResults',...
        'The gen_label is not consistent across mvpa results for averaging.');
elseif numel(unique({infos.gen_time})) > 1
    error('mvpa:averageMvpaRes:checkInfos:inconsistentMvpaResults',...
        'The gen_time is not consistent across mvpa results for averaging.');
elseif numel(unique({infos.sc_method})) > 1
    error('mvpa:averageMvpaRes:checkInfos:inconsistentMvpaResults',...
        'The sc_method is not consistent across mvpa results for averaging.');
elseif numel(unique({infos.svm_type})) > 1
    error('mvpa:averageMvpaRes:checkInfos:inconsistentMvpaResults',...
        'The svm_type is not consistent across mvpa results for averaging.');
elseif numel(unique({infos.tr_cond})) > 1
    error('mvpa:averageMvpaRes:checkInfos:inconsistentMvpaResults',...
        'The tr_cond is not consistent across mvpa results for averaging.');
elseif numel(unique({infos.tr_label})) > 1
    error('mvpa:averageMvpaRes:checkInfos:inconsistentMvpaResults',...
        'The tr_label is not consistent across mvpa results for averaging.');
elseif numel(unique({infos.tr_method})) > 1
    error('mvpa:averageMvpaRes:checkInfos:inconsistentMvpaResults',...
        'The tr_method is not consistent across mvpa results for averaging.');
elseif size(unique(cell2mat(cat(1,infos.tr_timePoints)),'rows'),1) > 1
    error('mvpa:averageMvpaRes:checkInfos:inconsistentMvpaResults',...
        'The tr_timePoints are not consistent across mvpa results for averaging.');
end

end

function avgData = averageIndividualData(objList,dataType,avgType)
% Function for averaging individual level data

temp = objList{1}.getInfo;
genTime = temp.gen_time;

switch dataType
    case 'gen_neuroMetrFun'
        if strcmp(avgType,'acrossSub')
            indivData = cellfun(@getNeuroMetricFunctions,objList,...
                                repmat({'genTime'},size(objList)),repmat({genTime},size(objList)));
            fieldExclusionStr = '';
        else
            predLabels = cellfun(@getPredLabels,objList,...
                                 repmat({'genTime'},size(objList)),repmat({genTime},size(objList)),...
                                 'UniformOutput',false);
            predLabelInfo = cellfun(@getInfoGenExamples,objList,...
                                    'UniformOutput',false);
            predLabels = cat(1,predLabels{:});
            predLabelInfo = cat(1,predLabelInfo{:});
            indivData = mvpa.fitNeurometricFunctions(predLabels,predLabelInfo,genTime);
            fieldExclusionStr = '';
        end
    case 'gen_neuroMetrFun_ffx'
        predLabels = cellfun(@getPredLabels,objList,...
            repmat({'genTime'},size(objList)),repmat({genTime},size(objList)),...
                            'UniformOutput',false);
        predLabelInfo = cellfun(@getInfoGenExamples,objList,...
                            'UniformOutput',false);
        predLabels = cat(1,predLabels{:});
        predLabelInfo = cat(1,predLabelInfo{:});
        indivData = mvpa.fitNeurometricFunctions(predLabels,predLabelInfo,genTime);
        fieldExclusionStr = '';
    case 'gen_perfEstimates'
        indivData = cellfun(@getGenPerfEstimates,objList,...
            repmat({'genTime'},size(objList)),repmat({genTime},size(objList)));
        fieldExclusionStr = '_perFolds$';
    case 'gen_recalIndex'
        indivData = cellfun(@getRecalIndex,objList,...
            repmat({'genTime'},size(objList)),repmat({genTime},size(objList)));
        fieldExclusionStr = '_err';
    case 'gen_recalIndex_ffx'
        predLabels = cellfun(@getPredLabels,objList,...
            repmat({'genTime'},size(objList)),repmat({genTime},size(objList)),...
                            'UniformOutput',false);
        predLabelInfo = cellfun(@getInfoGenExamples,objList,...
                            'UniformOutput',false);
        predLabels = cat(1,predLabels{:});
        predLabelInfo = cat(1,predLabelInfo{:});
        [indivData.RI,indivData.RI_err] = ...
            mvpa.compRecalIndex(predLabelInfo,predLabels);
        fieldExclusionStr = '';
end

avgData = indivData(1);
fieldNames = fieldnames(avgData);
fieldsToAverage = fieldNames(cellfun(@isempty,regexp(fieldNames,fieldExclusionStr)));
avgData = rmfield(avgData,fieldNames(...
    ~cellfun(@isempty,regexp(fieldNames,fieldExclusionStr))));

for i = 1:numel(fieldsToAverage)
    
    if ismember(fieldsToAverage{i},{'PF','stimLevels'})
        avgData.(fieldsToAverage{i}) = indivData(1).(fieldsToAverage{i});
    elseif strcmp(dataType,'gen_neuroMetrFun') && strcmp(avgType,'withinSub')
        break;
    elseif ~isempty(regexp(dataType,'_ffx$','once'))
        break;
    else
        
        dimToCat = find(size(indivData(1).(fieldsToAverage{i})) > 1,1,'last')+1;
        temp = cat(dimToCat,indivData.(fieldsToAverage{i}));
        
        if ~isempty(regexp(fieldsToAverage{i},'^RI','once'))
            avgData.(fieldsToAverage{i}) = mean(temp,dimToCat);
            if strcmp(avgType,'acrossSub')
                % Adding group level error estimates
                fieldName = [fieldsToAverage{i},'_std'];
                avgData.(fieldName) = std(temp,0,dimToCat);
                fieldName = [fieldsToAverage{i},'_err'];
                avgData.(fieldName) = (std(temp,0,dimToCat))./ ...
                    sqrt(size(temp,3));
            end
        elseif ~isempty(regexp(fieldsToAverage{i},'^r2_','once'))
            % Fisher transform first, then average when averaging correlation
            % coefficients.
            avgData.(fieldsToAverage{i}) = ...
                squeeze((tanh(mean(atanh(sqrt(temp)),dimToCat))).^2);
        elseif ~isempty(regexp(fieldsToAverage{i},'^r_','once'))
            % Fisher transform first, then average when averaging correlation
            % coefficients.
            avgData.(fieldsToAverage{i}) = ...
                squeeze(tanh(mean(atanh(temp),dimToCat)));
        elseif ~isempty(regexp(fieldsToAverage{i},'^(b|int|rf)_','once'))
            avgData.(fieldsToAverage{i}) = mean(temp,dimToCat);
            if strcmp(avgType,'acrossSub')
                % Adding group level error estimates
                fieldName = [fieldsToAverage{i},'_std'];
                avgData.(fieldName) = std(temp,0,dimToCat);
                fieldName = [fieldsToAverage{i},'_err'];
                avgData.(fieldName) = (std(temp,0,dimToCat))./ ...
                    sqrt(size(temp,3));
            end
        elseif ~isempty(regexp(fieldsToAverage{i},'(_pv|_pctr)$','once'))
            if strcmp(avgType,'acrossSub')
                avgData.(fieldsToAverage{i}) = temp;
            else
                avgData.(fieldsToAverage{i}) = mean(temp,dimToCat);
            end
        end
    end
end
avgData = orderfields(avgData);

end