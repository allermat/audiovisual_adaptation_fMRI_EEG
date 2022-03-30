function out = getData(obj)
% Method for accessing the data property
% 
% USAGE:
%   out = getInfo(obj)
% INPUT:
%   obj (object): mvpares object
% OUTPUT:
%   out (struct): the mvpa result dataset as a structure

% Copyright(C) 2017, Mate Aller
% allermat@gmail.com

if isa(obj.data,'matlab.io.MatFile')
    out = load(obj.source);
else
    out = obj.data;
end
    
end