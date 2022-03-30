function structArray = catpadstruct(varargin)
% Concatenates structures into a structure array, missing filelds padded with empy arrays

assert(all(cellfun(@isstruct,varargin)));

structArray = repmat(struct(),size(varargin));

for i = 1:numel(structArray)
    fn = fieldnames(varargin{i});
    for iField = 1:numel(fn)
        structArray(i).(fn{iField}) = varargin{i}.(fn{iField});
    end
end

end