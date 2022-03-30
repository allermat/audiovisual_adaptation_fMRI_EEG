function out = run_bold_fMRI(varargin)

% Parsing input
p = inputParser;

validHemispheres = {'mean_lhrh','lh','rh'};

addParameter(p,'hemisphere','mean_lhrh',@(x) ismember(x,validHemispheres));

parse(p,varargin{:});

hemisphere = p.Results.hemisphere;

% Get data
[gimg,aloc,conds,roiNames] = lme_getData('hemisphere',hemisphere);

% Assign data to output structure
out = struct;
out.gimg = gimg;
out.aloc = aloc;
out.conds = conds;
out.roiNames = roiNames;
