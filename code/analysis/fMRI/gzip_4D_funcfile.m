function gzip_4D_funcfile(funcfolder, expression)
%   Gzip 4D functional file

funcfile = cellfun(@(x) fullfile(x, getfname(x, expression)), funcfolder);
cellfun(@(x) gzip(x), funcfile);
delete(funcfile{:});