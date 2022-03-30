function Gc = getGmatrix(A,theta,nConds,varargin)

p = inputParser;

addRequired(p,'A',@(x) validateattributes(x,{'numeric'},{'finite'}));
addRequired(p,'theta',@(x) validateattributes(x,{'numeric'},{'finite'}));
addRequired(p,'nConds',@(x) validateattributes(x,{'numeric'},{'scalar','positive'}));
addParameter(p,'center',true,@(x) validateattributes(x,{'logical'},{'scalar'}));

parse(p,A,theta,nConds,varargin{:});

A = p.Results.A;
theta = p.Results.theta;
nConds = p.Results.nConds;
center = p.Results.center;

% Sum features weighted by thetas and center to get G matrix
M = zeros(nConds,size(A,2));
H = eye(nConds)-ones(nConds)/nConds;
for i = 1:size(A,3)
    M = M + theta(i)*A(:,:,i);
end
G = M*M';
if center
    Gc = H*G*H';
else
    Gc = G;
end
end