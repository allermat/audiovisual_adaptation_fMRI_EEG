function M = mdl2mat(mdl)
% Converts the libsvm's svmtrain output model to a matrix. 
% 
% DETAILS:
%   CAUTION! 
%       - Only tested with libsvm-3.20 
%       - Doesn't work with multi class classification. 
%       - Probability estimates are not saved. 
% 
% USAGE: 
%   M = mdl2mat(mdl)
%
% INPUT:
%   mdl: libsvm's svmtrain output structure. 
% 
% OUTPUT: 
%   M: the data of mdl converted to a matrix
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
%% Parse input
p = inputParser;

fieldsReqMdl = {'Parameters','nr_class','totalSV','rho','Label','sv_indices',...
    'ProbA','ProbB','nSV','sv_coef','SVs'};
checkMdl = @(x) all(isfield(x,fieldsReqMdl));

addRequired(p,'mdl',checkMdl);
parse(p,mdl);

mdl = p.Results.mdl;

%%

if mdl.nr_class > 2
    error('Multi-class classification is not supported!');
end

% If mdl.Parameters(1) > 1 mdl.Label and mdl.nSV is going to be empty
temp = [mdl.Parameters',mdl.nr_class,mdl.totalSV,mdl.rho,mdl.Label',mdl.nSV'];
% Matrix for storing mdl.Parameters, mdl.nr_class, mdl.totalSV, mdl.rho, 
% mdl.Label, mdl.nSV, mdl.sv_indices and mdl.sv_coef
M1 = NaN(3,max([mdl.totalSV,size(temp,2)]));
M1(1,1:size(temp,2)) = temp;
M1(2,1:size(mdl.sv_indices,1)) = mdl.sv_indices';
M1(3,1:size(mdl.sv_coef,1)) = mdl.sv_coef';
% Matrix for storing mdl.SVs
M2 = NaN(size(mdl.SVs,2),max([mdl.totalSV,size(temp,2)]));
M2(:,1:mdl.totalSV) = full(mdl.SVs');
% Concatenating the two matrices
M = [M1;M2];

end

