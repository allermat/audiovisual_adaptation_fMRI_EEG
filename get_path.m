function outpath = get_path(pathtype)

switch pathtype
    case {'main' 'root'}
        [outpath, dirname, ext] = fileparts(strrep(which(mfilename), [filesep mfilename '.m'], ''));
        
    case 'project'
        outpath = fileparts(mfilename('fullpath'));

    case 'misc'
        [outpath, fname, ext] = fileparts(which(mfilename));
        
    case 'toolbox'
        repo_root = fileparts(mfilename('fullpath'));
        outpath = fullfile(repo_root,'code','toolbox');
end
end