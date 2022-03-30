function sigma = compCovEEG(y,X)
% Method for computing covariance matrix on EEG data
% 
% USAGE:
%   sigma = compCovEEG(y,X)
% INPUT:
%   y (numeric): 3D array of data, format should be 
%       nExamples x nChannels x nTimePoints
%   X (numeric): Design matrix nExamples x nRegressors
% OUTPUT: 
%   sigma (numeric): 3D array of residuals, same size as y
% 

% Copyright(C) 2017, Mate Aller
% allermat@gmail.com

% Computing residuals
r = mvpa.compResidualsEEG(y,X);
df = size(X,1) - size(X,2);

sigma = NaN(size(y,2),size(y,2),size(y,3));

parfor iTimePoint = 1:size(y,3)
    sigma(:,:,iTimePoint) = rsa.stat.covdiag(r(:,:,iTimePoint),df);
end

end

