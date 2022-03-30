function r = compResidualsEEG(y,X)
% Method for computing residuals on EEG data based on least squares estimate
% 
% USAGE:
%   r = compResidualsEEG(y,X)
% INPUT:
%   y (numeric): 3D array of data, format should be 
%       nExamples x nChannels x nTimePoints
%   X (numeric): Design matrix nExamples x nRegressors
% OUTPUT: 
%   r (numeric): 3D array of residuals, same size as y
% 

% Copyright(C) 2017, Mate Aller
% allermat@gmail.com

r = NaN(size(y));
for iChannel = 1:size(y,2)
    parfor iTimePoint = 1:size(y,3)
        [~,~,r(:,iChannel,iTimePoint)] = ...
            regress(y(:,iChannel,iTimePoint),X);
    end
end

end

