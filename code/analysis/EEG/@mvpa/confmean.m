function [ciLow,ciHigh] = confmean(A,dim,alpha)
% Method for computing confidence interval of mean
% Standard error of the mean
SEM = std(A,0,dim)/sqrt(size(A,dim));
% T-Score
ts = tinv(1-alpha,size(A,dim)-1);
% Confidence Intervals 
ciLow = mean(A,dim) - ts*SEM;
ciHigh = mean(A,dim) + ts*SEM;
end

