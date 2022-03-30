function newPath = repairPath(oldPath)
% Repairs a saved path to a file so that it is valid in the current environment

% Parse input
% check if input is a valid path
checkPath = @(x) ~isempty(regexp(x,'([a-zA-Z_0-9:]+[\\/])([a-zA-Z_0-9.]+[\\/]?)*','once'));
p = inputParser;
addRequired(p,'oldPath',checkPath);
parse(p,oldPath);
oldPath = p.Results.oldPath;

currProjPath = strsplit(get_path('project'),filesep);

oldPathParts = strsplit(oldPath,'\');
idx = find(ismember(oldPathParts,currProjPath{end}),1);

if isempty(idx)
    error('The provided file path is not compatible with the current project');
end
newPathParts =  oldPathParts;
newPathParts(1:idx-1) = currProjPath(1:idx-1);

newPath = fullfile(newPathParts{:});

end