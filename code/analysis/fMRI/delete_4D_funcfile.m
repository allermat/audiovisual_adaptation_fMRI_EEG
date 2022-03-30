function delete_4D_funcfile(funcfolder, expression)
%   Delete 4D functional file that is not to be used any more

funcfile = cellfun(@(x) fullfile(x, getfname(x, expression)), funcfolder);
delete(funcfile{:});