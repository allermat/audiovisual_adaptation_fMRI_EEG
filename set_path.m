function set_path

% Get location of this file and use it as a base folder fo the repository
repo_root = fileparts(mfilename('fullpath'));
addpath(repo_root);

% Add toolbox folder to path
addpath(fullfile(repo_root,'code','toolbox'));

% All folders in toolbox to path
d = dir(fullfile(repo_root,'code','toolbox'));
dirNames = {d.name}';
dirNames = dirNames([d.isdir]);
dirNames = dirNames(~ismember(dirNames,{'.','..'}));
dirPaths = fullfile(repo_root,'code','toolbox',dirNames);
if isempty(dirPaths)
    warning('Toolboxes might be missing. Please refer to the dependencies in README.')
else
    addpath(dirPaths{:});
    % Initialize fieldtrip if it is added to the path
    if ismember('fieldtrip',dirNames)
        ft_defaults;
    end
    
    % libsvm's matlab subfolder must be added to the path
    m = regexp(dirNames,'libsvm-[0-9.]+','match','once');
    if any(~cellfun(@isempty,m))
        libsvm_dir = m{~cellfun(@isempty,m)};
        addpath(fullfile(repo_root,'code','toolbox',libsvm_dir,'matlab'));
    end
end

% Add analysis folders
addpath(fullfile(repo_root,'code','analysis','psychophysics'));
addpath(fullfile(repo_root,'code','analysis','fMRI'));
addpath(fullfile(repo_root,'code','analysis','EEG'));
% addpath(fullfile(path_repo, 'eyelink'));
addpath(fullfile(repo_root,'results'));
addpath(fullfile(repo_root,'results','data'));

end