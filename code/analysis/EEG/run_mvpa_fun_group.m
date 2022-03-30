function varargout = run_mvpa_fun_group(fun,trMethod,matchStr,varargin)
% Runs mvpa functions on the group of subjects

%% Parsing input, checking matlab
p = inputParser;

validTrMethods = {'sample-wise','sample-wise-avg'};

addRequired(p,'fun',@(x) validateattributes(x,{'function_handle'},{'scalar'}));
addRequired(p,'trMethod',@(x) any(validatestring(x, ...
                                                 validTrMethods)));
addRequired(p,'matchStr',@ischar);
addParameter(p,'funInput',{},@(x) validateattributes(x,{'cell'},{'vector','row'}));
addParameter(p,'nFilesExpected',[],@(x) validateattributes(x,{'numeric'},{'scalar'}));
addParameter(p,'subFolder',[],@(x) validateattributes(x,{'char'},{'nonempty'}));

parse(p,fun,trMethod,matchStr,varargin{:});

fun = p.Results.fun;
trMethod = p.Results.trMethod;
fileMatchStr = p.Results.matchStr;
funInput = p.Results.funInput;
nFilesExpected = p.Results.nFilesExpected;
subFolder = p.Results.subFolder;

% Finding subject folders
saveDf = cd(fullfile(get_path('project'),'data'));
subjList = cellstr(ls);
subjList = regexp(subjList,'[0-9]{3}','match','once');
subjList = subjList(~ismember(subjList,{''}));
cd(saveDf);

% Collecting subject level mvpares datasets
miscInput = {};
% Preparations specific to certain functions
if strcmp(func2str(fun),'mvpa.mergeMvpaRes')
    % Expected Number of files per subject
    if isempty(nFilesExpected)
        nFilesExpected = 2;
    end
    filePathList = collectFiles(subjList,fileMatchStr,trMethod,nFilesExpected,subFolder);
else
    % Expected Number of files per subject
    if isempty(nFilesExpected)
        nFilesExpected = 1;
    end
    filePathList = collectFiles(subjList,fileMatchStr,trMethod,nFilesExpected,subFolder);
end

% Executing function on all subject level data
if strcmp(func2str(fun),'mvpa.averageMvpaRes')
    id = regexp(filePathList{1},'.*_(tr[A-Za-z-() ]+_).*(gen[A-Za-z-() ]+)', ...
                'tokens');
    id = id{:};
    id = cat(2,id{:});
    % Getting rid of special caracters to create valid file name
    id = regexprep(id, '[/\*:?"<>|]*',' ');
    avgFileName = ['gr_',id,'.mat'];
    
    I.pathAveragedFile = fullfile(get_path('project'), ...
                                  'data','group','EEG','MVPA',trMethod,avgFileName);
    if ~exist(fullfile(get_path('project'),'data','group','EEG','MVPA',trMethod),'file')
        mkdir(fullfile(get_path('project'),'data','group','EEG','MVPA',trMethod));
    end
    I.pathFilesToAverage = filePathList;
    varargout{1} = mvpa.averageMvpaRes(I);
elseif strcmp(func2str(fun),'mvpa.mergeMvpaRes')
    id = regexp(fileMatchStr,'([cenr]+_)(tr[A-Za-z-()|?]+_).*_(gen[A-Za-z-()|?]+)', ...
                'tokens');
    id = id{:};
    id = cat(2,id{:});
    % Getting rid of special caracters to create valid file name
    id = regexprep(id, '[/\*:?"<>|]*',' ');
    for i = 1:numel(subjList)
        
        subID = subjList{i};
        
        mergedFileName = [id,'.mat'];
        
        I.pathMergedFile = fullfile(get_path('project'), ...
                                  'data',subID,'EEG','MVPA',trMethod,mergedFileName);
        I.pathFilesToMerge = filePathList(i,:);
        varargout{i} = mvpa.mergeMvpaRes(I); %#ok
        
    end
elseif strcmp(func2str(fun),'collect')
    for i = 1:size(filePathList,1)
        % try 
            % temp{i} = mvpares(filePathList{i}); %#ok
        % catch
            temp(i) = load(filePathList{i}); %#ok
        % end
    end
    varargout{1} = temp;
else
    for i = 1:size(filePathList,1)
        temp = mvpares(filePathList{i});
        if ~isempty(miscInput)
            % This cell array must be a row vector
            args = cat(2,funInput,miscInput(i,:));
        else
            args = funInput;
        end
        varargout{i} = fun(temp,args{:}); %#ok
    end
end

end

function filePathList = collectFiles(subjList,fileMatchStr,trMethod,nFilesExpected,subFolder)

filePathList = {};
index = 1;
for i = 1:size(subjList,1)
    saveDf = cd(fullfile(get_path('project'),'data',subjList{i},'EEG','MVPA',trMethod,subFolder));
    fileList = cellstr(ls);
    matchID = ~cellfun(@isempty,regexp(fileList,fileMatchStr));
    if sum(matchID) == 0
        warning('No file, skipping subject %s! ',subjList{i});
        cd(saveDf);
        continue;
    elseif sum(matchID) > nFilesExpected
        warning('More files than needed, skipping subject %s! ',subjList{i});
        cd(saveDf);
        continue;
    else
        fileName = fileList(matchID);
        if iscolumn(fileName)
            fileName = fileName';
        end
    end
    temp = cellfun(@fullfile,repmat({pwd},size(fileName)),fileName,'UniformOutput',false);
    filePathList = cat(1,filePathList,temp);
    index = index + 1;
    cd(saveDf);
end

if isrow(filePathList)
    filePathList = filePathList';
end

end
