function [y,RDM_euc_vox] = decision_uncertainty_model(varargin)
% Function for the decision model
% 
% USAGE:
%   population_rate_code_model('Name',value)
% INPUT:
%   'Name',value pair arguments:
%       aloc: vector of auditory locations of interest
%       nNeurons: number of simulated neurons/voxels
%           default = 500;
%       propResp: proportion of neurons responding to uncertainty
%           default = 0.7
%       sigma: standard deviation of underlying PF in degrees
%           default = 10
%       mu: location of maximum uncertainty in degrees, default = 90
%       span: spread of the sigmas across voxels in
%           degrees, default = 20
%       shift: lateral shift in tuning curves in degrees, default = 0
%       plotFigure: whether or not to plot output figure
%       seed: seed for the random number generator
%   
% OUTPUT:
%   y: simulated neuron level responses 
%       (nNeurons x nConditions)
%   Optional figure with raw and sampled data images and RDMs
%

% Copyright (C) 2019, Mate Aller
% allermat@gmail.com

% Parsing input
p = inputParser;
addParameter(p,'aloc',[-12,-5,-2,0,2,5,12],@(x) validateattributes(x,{'numeric'},...
                                                  {'row','finite'}));
addParameter(p,'nNeurons',1e5,@(x) validateattributes(x,{'numeric'},...
                                                  {'scalar','positive','integer'}));
addParameter(p,'propResp',0.7,@(x) validateattributes(x,{'numeric'},...
                                                  {'scalar','>=',0,'<=',1}));
addParameter(p,'sigma',10,@(x) validateattributes(x,{'numeric'}, ...
                                                {'scalar','nonnegative','finite'}));
addParameter(p,'mu',0,@(x) validateattributes(x,{'numeric'}, ...
                                                {'scalar','nonnegative','<=',360}));
addParameter(p,'span',5,@(x) validateattributes(x,{'numeric'}, ...
                                                {'scalar','finite','nonnegative'}));
addParameter(p,'shift',0,@(x) validateattributes(x,{'numeric'}, ...
                                                 {'scalar','finite'}));
addParameter(p,'seed',[],@(x) validateattributes(x,{'numeric'},...
                                                  {'scalar','finite','nonnegative'}));
addParameter(p,'plotFigure',true,@(x) validateattributes(x,{'logical'},{'scalar'}));

parse(p,varargin{:})

aloc = p.Results.aloc;
nNeurons = p.Results.nNeurons;
propResp = p.Results.propResp;
sigma = p.Results.sigma;
mu = p.Results.mu;
span = p.Results.span;
shift = p.Results.shift;
plotFigure = p.Results.plotFigure;
seed = p.Results.seed;

% Further specific input check
if sigma < span/2
    error('decision_model:invalidInput',...
          'Span to high, this would result in negative sigma');
end

% Reseeding the random number generator
if isempty(seed)
    rng('shuffle');
else
    rng(seed);
end

% Uncertainty function converted from PF using eq. 2 from 
% Grinband et al. 2006 except instead of a Weibull function I use a
% cumulative gaussian. This has a maximum of 1 at the PSE of the PF and 
% diminishes towards the ends
fun_UCF = @(x) abs(abs(x-0.5)-0.5)*2;
% Voxels responding to uncertainty
y = fun_UCF(normcdf(repmat(aloc+shift,round(nNeurons*propResp),1),mu,...
            repmat(sigma+(rand(round(nNeurons*propResp),1)*span)-(span/2),1,size(aloc,2))));
% Voxels not responding to uncertainty
temp = rand(round(nNeurons*(1-propResp)),numel(aloc));
y = cat(1,y,temp); 

% Computing RDM for raw data
RDM_euc_vox = squareform(pdist(y','euclidean'));

% % Computing RDM for sampled data
% RDM_euc_sampl = squareform(pdist(y_sampl','squaredeuclidean'));

% Plotting
if plotFigure
    figure;
    set(gcf, 'Position', [100 600 850 300]);

    subplot(1,2,1);
    imagesc(y);
    line([0,numel(aloc)+1],[nNeurons(1),nNeurons(1)],'Color','w','LineWidth',2);
    set(gca,'Xtick',1:numel(aloc),'XTickLabels',aloc,'Ytick',[]);
    title('Voxel responses, samples x locations');

    subplot(1,2,2);
    image(rsa.util.scale01(rsa.util.rankTransform_equalsStayEqual(RDM_euc_vox,1)), ...
          'CDataMapping','scaled');
    set(gca,'CLim',[0 1],'CLimMode','manual');
    set(gca,'Xtick',1:numel(aloc),'XTickLabels',num2cell(aloc),...
            'Ytick',1:numel(aloc),'YTickLabels',num2cell(aloc));
    title('RDM of voxel responses (Euclidean)');
end

end