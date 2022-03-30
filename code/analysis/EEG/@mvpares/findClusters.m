function clusters = findClusters(h)
% Method for finding clusters of true values in input vector
% 
% USAGE:
%   clusters = findClusters(h,time)
% INPUT: 
%   h (logical): vector of logicals
% OUTPUT:
%   clusters (numerical): N x 2 matrix of start and end time values of
%       clusters. N = number of clusters, column 1: start indices, 
%       column 2: end indices. If there are no clusters it returns and
%       empty array.

% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

% Parsing input
p = inputParser;
addRequired(p,'h',@(x) validateattributes(x,{'logical'},{'vector'}));
parse(p,h);
h = p.Results.h;

if isrow(h), h = h'; end

hIdx = find(h);
if any(hIdx)
    temp = diff(hIdx);
    temp = [2;temp];
    startIdx = find(temp > 1);
    endIdx = [startIdx(2:end)-1;size(temp,1)];
    clusters = [hIdx(startIdx),hIdx(endIdx)];
else
    clusters = [];
end

end

