function funcfile = get_3D_funcfile(funcfolder, expression, subject)
%   Gunzip 4D functional file and add 3D volume indexes

% Number of folders
nfolders = size(funcfolder, 1);

% Extra volumes to be deleted
extravol = zeros(nfolders, 1);
if exist('subject', 'var') && isfield(subject.fMRI, 'extravol')
    extravolfile = fullfile(get_path('project'), get_folder(subject, 'r', 'fMRI', 'processed data', 'extravol'));
    if size(extravolfile, 1) == 1
        id = ~cellfun(@isempty, strfind(funcfolder, extravolfile));
    else
        id = any(cell2mat(cellfun(@(x) ~cellfun(@isempty, strfind(funcfolder, x)), extravolfile', 'un', 0)), 2);
    end
    extravol(id) = cat(1, subject.fMRI.extravol.delete{:});
end

funcfile = cellfun(@(x) fullfile(x, getfname(x, expression)), funcfolder);
for s=1:nfolders
    if strfind(expression, 'gz')
        funcfile(s) = gunzip(funcfile{s}); % get unzipped 4D functional
    end
    h = spm_vol(funcfile{s});
    nvols = numel(h) - extravol(s); % delete extra volumes if needed
    funcfile{s} = cellfun(@(x,y) [x sprintf(',%d', y)], repmat(funcfile(s), nvols, 1), num2cell((1:nvols))', 'un', 0); % append volume index/frame
end
