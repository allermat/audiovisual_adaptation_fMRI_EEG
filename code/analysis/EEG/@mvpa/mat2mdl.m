function mdl = mat2mdl(M)
% Converts the libsvm's svmtrain output model converted to a matrix back to sturcture form. 
% 
% DETAILS:
%   CAUTION! 
%       - Only tested with libsvm-3.20 
%       - Doesn't work with multi class classification. 
%       - Probability estimates are not saved.
% 
% USAGE: 
%   mdl = mat2mdl(M)
%
% INPUT:
%   M: the data of libsvm's svmtrain output structure converted to a matrix
%       Dimensions of the matrix: 
%           nRows: nFeatures+3 
%           nCols: number of totalSV or 8 or 12 depending on SV type and 
%                  number of totalSV
%       M(1,1:5) = mdl.Parameters'
%       M(1,6) = mdl.nr_class
%       M(1,7) = mdl.totalSV
%       M(1,8) = mdl.rho
%       M(1,9:10) = mdl.Label' - if 2 class classification, otherwise NaN
%       M(1,11:12) = mdl.nSV' - if 2 class classification, otherwise NaN
%       M(2,:) = mdl.sv_indices'
%       M(3,:) = mdl.sv_coef'
%       M(4:end,:) = mdl.SVs'
% 
% OUTPUT: 
%   mdl: libsvm's svmtrain output structure. 
% 
%% Parse input
p = inputParser;

addRequired(p,'mdl',@ismatrix);
parse(p,M);

M = p.Results.mdl;

%%

if M(1,6) > 2
    error('Multi-class classification is not supported!');
end

% Generating model structure
mdl = struct();
mdl.Parameters = M(1,1:5)';
mdl.nr_class = M(1,6);
mdl.totalSV = M(1,7);
mdl.rho = M(1,8);
if M(1,1) > 1
    mdl.Label = [];
else
    mdl.Label = M(1,9:10)';
end
mdl.sv_indices = M(2,1:mdl.totalSV)';
mdl.ProbA = [];
mdl.ProbB = [];
if M(1,1) > 1
    mdl.nSV = [];
else
    mdl.nSV = M(1,11:12)';
end
mdl.sv_coef = M(3,1:mdl.totalSV)';
mdl.SVs = sparse(M(4:end,1:mdl.totalSV)');


end

