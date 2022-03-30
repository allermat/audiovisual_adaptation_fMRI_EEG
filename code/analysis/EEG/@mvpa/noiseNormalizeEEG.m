function y = noiseNormalizeEEG(y,s)
% Method for noise normalization on EEG data
% 
% USAGE:
%   sigma = compCovEEG(y,X)
% INPUT:
%   y (numeric): 3D array of data, format should be 
%       nExamples x nChannels x nTimePoints
%   s (numeric): covariance matrices, format should be 
%       nChannels x nChannels x nTimePoints
% OUTPUT: 
%   y (numeric): 3D array of activity patterns after noise
%       normalization
% 

% Copyright(C) 2017, Mate Aller
% allermat@gmail.com
% The code below uses a section from the rsatoolbox function
% noiseNormalizeBeta written by
% Alexander Walther, Joern Diedrichsen
% joern.diedrichsen@googlemail.com

if size(s,3) == 1
    [V,L]=eig(s);
    l=diag(L);
    sq = V*bsxfun(@rdivide,V',sqrt(l)); % Slightly faster than sq = V*diag(1./sqrt(l))*V';
    parfor i = 1:size(y,3)
        % Postmultiply by the inverse square root of the estimated
        % matrix. Code snippet from the rsatoolbax
        y(:,:,i)=y(:,:,i)*sq;
    end
else
    parfor i = 1:size(y,3)
        % Postmultiply by the inverse square root of the estimated
        % matrix. Code snippet from the rsatoolbax
        [V,L]=eig(s(:,:,i));
        l=diag(L);
        sq = V*bsxfun(@rdivide,V',sqrt(l)); % Slightly faster than sq = V*diag(1./sqrt(l))*V';
        y(:,:,i)=y(:,:,i)*sq;
    end
end
end

