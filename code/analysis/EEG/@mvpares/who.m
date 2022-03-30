function out = who(obj)
% Method to keep the functionality of the Matlab.io.matFile object's who method
% 
% USAGE:
%   out = who(obj)
% INPUT:
%   obj (object): mvpares object
% OUTPUT: 
%   out (cell array): list of fields of the dataset

% Copyright(C) 2017, Mate Aller
% allermat@gmail.com

if isa(obj.data,'matlab.io.MatFile')
    out = who(obj.data);
else
    out = fieldnames(obj.data);
end

end