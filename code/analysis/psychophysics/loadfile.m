function dataout = loadfile(folder, substr, fileid, var)
%  dataout = loadfile(folder, substr, fileid, datain)
%  loadfile  loads and concatenates dataset(s) from file(s)

% Check the dimension of input parameters
nfolders = length(folder);
if nfolders ~= length(substr) || nfolders ~= length(fileid)
   error('Dimension of input parameters are not identical') 
end

% Put requested files together
for i=1:nfolders
    if verLessThan('matlab', '8.1')
        filename = getfname(folder{i}, substr{i});
        for j=1:length(filename)
            fname{j} = fullfile(folder{i}, filename{j});
        end
    else
        fname = fullfile(folder{i}, getfname(folder{i}, substr{i}));
    end
    if isempty(fname)
        error('The requested file is not found');
    end
    for j=1:length(fileid{i})
        files(j,i) = fname(fileid{i}(j));
    end
end
files = files(:); % concatenate all files
files = files(~cellfun(@isempty, files)); % remove empty cells
nfiles = length(files);

% Concatenate datasets
dataout = dataset;
for i=1:nfiles
    load(files{i}, var{:});
    dataout = [dataout; data];
end


