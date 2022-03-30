function [onset, names, blockstart] = load_onset(subject, regressor_unit)

behavfolder = fullfile(get_path('project'), get_folder(subject, 'r', 'fMRI', 'behavioural data'));
if ~iscell(behavfolder)
   behavfolder = {behavfolder}; 
end

for s=1:numel(behavfolder)
    fname{1} = arrayfun(@(x) getfname(behavfolder{s}, sprintf('*run%02d*', x)), subject.fMRI.pretest.runid{s}'); % pre-test filenames
    postdir = [lower(subject.fMRI.pretest.match{s}) 'adapt'];
    fname{2} = arrayfun(@(x) getfname(behavfolder{s}, sprintf('*run%02d*', x)), ...
        subject.fMRI.posttest.(postdir).runid{ismember(subject.fMRI.posttest.(postdir).session, subject.fMRI.pretest.session(s))}'); % corresponding post-test filenames
    fname = cat(1, fname{:}); % concatenate filenames
    for f=1:numel(fname)
       [data, time] = load_truncated_data(fullfile(behavfolder{s}, fname{f}), 20, 'data', 'time');
       T = dataset2table(data);
       T.onset = mean([time.vonset time.aonset], 2) - time.start; % add onset times
       T(isnan(T.aloc),:) = []; % remove fixation periods (must be after onset had been added!!!)
       blockstart{f,s} = T.onset(T.blockstart == 1);
       switch regressor_unit
           case 'block'
               blocks = unique(T.block);
               for b=1:numel(blocks)
                   blockdata = T(T.block == b,:);
                   blocktype = unique(blockdata.blocktype);
%                    blocktype{:}
                   if strfind(blocktype{:}, 'test')
                       aloc = unique(blockdata.aloc(~isnan(blockdata.aloc)));
                       nlocs = numel(aloc);
                       for i=1:nlocs
                           onset{f,s}{1,i,b} = blockdata.onset(blockdata.aloc == aloc(i) & isnan(blockdata.resp));
                           names{f,s}{1,i,b} = sprintf('session=%d,run=%d,block=%d,resp=0,aloc=%d', unique(blockdata.session(blockdata.block == blocks(b))), ...
                               unique(blockdata.run(blockdata.block == blocks(b))), blocks(b), aloc(i));
                           onset{f,s}{2,i,b} = blockdata.onset(blockdata.aloc == aloc(i) & ~isnan(blockdata.resp));
                           names{f,s}{2,i,b} = sprintf('session=%d,run=%d,block=%d,resp=1,aloc=%d', unique(blockdata.session(blockdata.block == blocks(b))), ...
                               unique(blockdata.run(blockdata.block == blocks(b))), blocks(b), aloc(i));
                       end
                   elseif strfind(blocktype{:}, 'adapt')
                       blockdata.v_blockstart = zeros(size(blockdata, 1), 1);
                       blockdata.v_blockstart(unique([1; find(blockdata.vloc ~= circshift(blockdata.vloc, 1))])) = 1; % different vloc-s are modelled!! (block starts with new vloc)
                       vloc = unique(blockdata.vloc(~isnan(blockdata.vloc)));
                       nlocs = numel(vloc);
                       for i=1:nlocs
                           onset{f,s}{1,i,b} = blockdata.onset(blockdata.v_blockstart & blockdata.vloc == vloc(i)); % & isnan(blockdata.resp)
                           names{f,s}{1,i,b} = sprintf('session=%d,run=%d,block=%d,resp=0,vloc=%d', unique(blockdata.session(blockdata.block == blocks(b))), ...
                               unique(blockdata.run(blockdata.block == blocks(b))), blocks(b), vloc(i));
                       end
                       onset{f,s}{1,nlocs+1,b} = blockdata.onset(ismember(blockdata.vloc, vloc) & ~isnan(blockdata.resp));
                       names{f,s}{1,nlocs+1,b} = sprintf('session=%d,run=%d,block=%d,resp=1,vloc', unique(blockdata.session(blockdata.block == blocks(b))), ... % all V responses modelled separately!!
                           unique(blockdata.run(blockdata.block == blocks(b))), blocks(b));
                   else
                       error('unknown blocktype found');
                   end
               end
           case 'run'
               testdata = T(~cellfun(@isempty, strfind(T.blocktype, 'test')),:);
               aloc = unique(testdata.aloc(~isnan(testdata.aloc)));
               nlocs = numel(aloc);
               for i=1:nlocs
                   onset{f,s}{1,i,1} = testdata.onset(testdata.aloc == aloc(i) & isnan(testdata.resp));
                   names{f,s}{1,i,1} = sprintf('session=%d,run=%d,resp=0,aloc=%d', unique(testdata.session), unique(testdata.run), aloc(i));
                   onset{f,s}{2,i,1} = testdata.onset(testdata.aloc == aloc(i) & ~isnan(testdata.resp));
                   names{f,s}{2,i,1} = sprintf('session=%d,run=%d,resp=1,aloc=%d', unique(testdata.session), unique(testdata.run), aloc(i));
               end
               adaptdata =  T(~cellfun(@isempty, strfind(T.blocktype, 'adaptation')),:);
               if ~isempty(adaptdata)
                   adaptdata.v_blockstart = zeros(size(adaptdata, 1), 1);
                   adaptdata.v_blockstart(unique([1; find(adaptdata.vloc ~= circshift(adaptdata.vloc, 1))])) = 1; % different vloc-s are modelled!! (block starts with new vloc)
                   vloc = unique(adaptdata.vloc(~isnan(adaptdata.vloc)));
                   nlocs = numel(vloc);
                   for i=1:nlocs
                       onset{f,s}{1,i,2} = adaptdata.onset(adaptdata.v_blockstart & adaptdata.vloc == vloc(i));
                       names{f,s}{1,i,2} = sprintf('session=%d,run=%d,resp=0,vloc=%d', unique(adaptdata.session), unique(adaptdata.run), vloc(i));
                   end
                   onset{f,s}{1,nlocs+1,2} = adaptdata.onset(ismember(adaptdata.vloc, vloc) & ~isnan(adaptdata.resp));
                   names{f,s}{1,nlocs+1,2} = sprintf('session=%d,run=%d,resp=1,vloc', unique(adaptdata.session), unique(adaptdata.run)); % all V responses modelled separately!!
               end
       end
       onset{f,s} = reshape(onset{f,s}, 1, []);
       names{f,s} = reshape(names{f,s}, 1, []);
       onset{f,s}(cellfun(@isempty, names{f,s})) = []; % remove not modelled blocks
       names{f,s}(  cellfun(@isempty, names{f,s})) = []; % remove not modelled blocks
    end
    clear fname;
end
onset = reshape(onset, [], 1);
onset(cellfun(@isempty, onset)) = [];
names = reshape(names, [], 1);
names(cellfun(@isempty, names)) = [];
blockstart = reshape(blockstart, [], 1);
blockstart(cellfun(@isempty, blockstart)) = [];
