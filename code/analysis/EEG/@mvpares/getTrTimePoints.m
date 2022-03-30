function trTimePoints = getTrTimePoints(obj)
% Method for getting the training time points
%
% Copyright(C) 2016, Mate Aller
if any(cellfun(@numel,obj.info.tr_timePoints) > 1)
    fun = str2func(obj.info.timeBinRegFun);
    trTimePoints = cellfun(fun,obj.info.tr_timePoints)';
else
    trTimePoints = cell2mat(obj.info.tr_timePoints)';
end
end