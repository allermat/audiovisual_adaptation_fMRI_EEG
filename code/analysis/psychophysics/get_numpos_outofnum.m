function [NumPos, OutOfNum] = get_numpos_outofnum(data)
% Calculate number of positive responses and out of number responses from
% data

[NumPos, OutOfNum] = deal(cell(size(data)));
for i=1:numel(data)
    if ~isempty(data{i})
        subdata = data{i}(~cellfun(@isempty, strfind(data{i}.blocktype, 'test')),:);
        StimLevels = unique(subdata.aloc(~isnan(subdata.aloc)))';
        for j=1:length(StimLevels)
            NumPos{i}(j) = sum(subdata.resp == 2 & subdata.aloc==StimLevels(j) & ~isnan(subdata.resp) & ~isnan(subdata.aloc));
            OutOfNum{i}(j) = sum(subdata.aloc==StimLevels(j) & ~isnan(subdata.resp) & ~isnan(subdata.aloc));
        end
    end
end

% Concatenate data
NumPos = cat(1, NumPos{:});
OutOfNum = cat(1, OutOfNum{:});